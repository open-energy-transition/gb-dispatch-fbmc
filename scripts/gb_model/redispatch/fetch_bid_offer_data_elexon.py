# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Fetch Bid Offer Data from Elexon BMRS API
"""

import logging
from pathlib import Path

import requests
import pandas as pd
from datetime import datetime
import numpy as np
import pypsa

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def fetch_api_request_data(url: str, retrieval_message: str) -> pd.DataFrame:
    """
    Fetch Bid Offer Data (BOD) from Elexon API for the specified date range.

    Parameters
    ----------
    url: str
        API request URL
    retrieval_message: str
        Message to log for the retrieval operation
    """

    df = pd.DataFrame()

    max_retries = 3
    for i in range(1,max_retries+1):

        try:
            headers = {'User-Agent': 'Mozilla/5.0'}
            
            response = requests.get(url, headers=headers)
            response.raise_for_status()

            json_data = response.json()

            # Convert JSON data to Dataframe
            df = pd.DataFrame(json_data) 
            
            logger.info(f"Successfully retrieved {retrieval_message}.")
            break # Exit retry loop if successful

        except requests.exceptions.RequestException as e:

            logger.warning(f" Attempt {i}/{max_retries}: Error in fetching {retrieval_message}: {e}")
            if i == max_retries:
                logger.error(f"Max retries reached. Failed to fetch {retrieval_message}.")
        
    return df


def get_year_bod(base_url: str, bmu_units_filter:str, year: int) -> pd.DataFrame:
    """
    Fetch bid / offer data for an entire year by querying day by day

    Parameters:
    ----------
    base_url: str
        Base URL for the Elexon BMRS API
    bmu_units_filter: str
        Filter string to specify BMU units in the API request
    year: int
        Year for which to fetch the data
    """
    dfs = []

    start = datetime(year, 1, 1)
    end = datetime(year, 12, 31)

    current = start

    while current <= end:
        # Retrieving data month by month
        next_month = current + pd.DateOffset(months=1)

        # API request URL for the month
        url = f"{base_url}/datasets/BOD/stream?from={current}&to={next_month}{bmu_units_filter}"

        # Fetch data for the current month
        df = fetch_api_request_data(url, retrieval_message=f"Bid/Offer data for the dates {current} to {next_month}")

        dfs.append(df)

        # Reset current date to next month
        current = next_month

    return pd.concat(dfs, ignore_index=True)


def fetch_BM_units(base_url: str, technology_mapping: dict[str, str]):
    """
    Fetch BM Unit data from Elexon API to get mapping of BM units to fuel types

    Parameters
    ----------
    base_url: str
        Base URL for the Elexon BMRS API
    """

    # URL to fetch BM unit data
    url = f"{base_url}/reference/bmunits/all"

    # Fetch BM unit data
    df_bmu = fetch_api_request_data(url, retrieval_message="BMU Unit Data")

    # Map BM unit fuel types to PyPSA carriers
    df_bmu['carrier'] = df_bmu['fuelType'].map(technology_mapping)

    # Filter to only include BM units relevant to conventional technologies in the model
    # Done to improve API request performance as CfD prices are used for renewables
    # Helps reduce memory requirements for storing the BOD data
    df_bmu_filtered = df_bmu.loc[df_bmu['carrier'].isin(technology_mapping.values())]
    bmu_units = df_bmu_filtered['nationalGridBmUnit'].unique().tolist()

    # Create filter string for API request to only include relevant BM units
    filter_bmu_units = ",".join([(f"&bmUnit={x}") for x in bmu_units]).replace(",","")

    logger.info(f"Created BM unit filter for API request with {len(bmu_units)} units.")
    return filter_bmu_units


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    base_url = "https://data.elexon.co.uk/bmrs/api/v1"

    filter_bmu_units = fetch_BM_units(base_url, technology_mapping=snakemake.params.technology_mapping)

    df_bod = get_year_bod(base_url, filter_bmu_units, year=int(snakemake.wildcards.bod_year))

    df_bod.to_csv(snakemake.output.csv, index=False)
    logger.info(f"Saved Bid/Offer data to {snakemake.output.csv}")