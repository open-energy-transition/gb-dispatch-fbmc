# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Fetch Bid Offer Data from Elexon BMRS API
"""

import logging
from datetime import datetime
from pathlib import Path

import pandas as pd
import requests

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
    for i in range(1, max_retries + 1):
        try:
            headers = {"User-Agent": "Mozilla/5.0"}

            response = requests.get(url, headers=headers)
            response.raise_for_status()

            json_data = response.json()

            # Convert JSON data to Dataframe
            df = pd.DataFrame(json_data)

            logger.info(f"Successfully retrieved {retrieval_message}.")
            break  # Exit retry loop if successful

        except requests.exceptions.RequestException as e:
            logger.warning(
                f" Attempt {i}/{max_retries}: Error in fetching {retrieval_message}: {e}"
            )
            if i == max_retries:
                logger.error(
                    f"Max retries reached. Failed to fetch {retrieval_message}."
                )

    return df


def get_year_bod(
    base_url: str, bmu_units_filter: str, bmu_carrier_map: dict[str, str], year: int
) -> pd.DataFrame:
    """
    Fetch bid / offer data for an entire year by querying day by day

    Parameters
    ----------
    base_url: str
        Base URL for the Elexon BMRS API
    bmu_units_filter: str
        Filter string to specify BMU units in the API request
    bmu_carrier_map: dict[str, str]
        Mapping of BM units to carriers
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
        df = fetch_api_request_data(
            url,
            retrieval_message=f"Bid/Offer price data for the dates {current} to {next_month}",
        )

        dfs.append(df)

        # Reset current date to next month
        current = next_month

    df_bod = pd.concat(dfs, ignore_index=True)

    df_bod["carrier"] = df_bod["nationalGridBmUnit"].map(bmu_carrier_map)

    df_bod_mean = df_bod.groupby("carrier")[["bid", "offer"]].mean()
    return df_bod_mean


def fetch_BM_units(
    base_url: str,
    technology_mapping: dict[str, str],
    bmu_fuel_map_path: str,
    api_bmu_fuel_map: bool,
) -> tuple[str, dict]:
    """
    Fetch BM Unit data from Elexon API to get mapping of BM units to fuel types

    Parameters
    ----------
    base_url: str
        Base URL for the Elexon BMRS API
    technology_mapping: dict[str, str]
        Map Elexon carrier types to PyPSA carrier types
    bmu_fuel_map_path: str
        CSV path for mapping of BMU units to their fueltype
    api_bmu_fuel_map: bool
        Boolean to choose between fetching BM units via API (True) / reading existing excel (False)
    """

    if api_bmu_fuel_map:
        # URL to fetch BM unit data
        url = f"{base_url}/reference/bmunits/all"

        # Fetch BM unit data
        df_bmu = fetch_api_request_data(url, retrieval_message="BMU Unit Data")
    else:
        df_bmu = pd.read_excel(bmu_fuel_map_path).rename(
            columns={"NESO BMU ID": "nationalGridBmUnit", "REG FUEL TYPE": "fuelType"}
        )

    # Map BM unit fuel types to PyPSA carriers
    df_bmu["carrier"] = df_bmu["fuelType"].map(technology_mapping)

    # Filter to only include BM units relevant to conventional technologies in the model
    # Done to improve API request performance as CfD prices are used for renewables
    # Helps reduce memory requirements for storing the BOD data
    df_bmu_filtered = df_bmu.loc[df_bmu["carrier"].isin(technology_mapping.values())]
    bmu_units = df_bmu_filtered["nationalGridBmUnit"].unique().tolist()

    # Create filter string for API request to only include relevant BM units
    filter_bmu_units = ",".join([(f"&bmUnit={x}") for x in bmu_units]).replace(",", "")

    # Dictionary to map BM units to carriers
    bmu_carrier_map = dict(df_bmu_filtered[["nationalGridBmUnit", "carrier"]].values)

    logger.info(f"Created BM unit filter for API request with {len(bmu_units)} units.")
    return filter_bmu_units, bmu_carrier_map


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    base_url = "https://data.elexon.co.uk/bmrs/api/v1"

    filter_bmu_units, bmu_carrier_map = fetch_BM_units(
        base_url,
        technology_mapping=snakemake.params.technology_mapping,
        bmu_fuel_map_path=snakemake.input.bmu_fuel_map_path,
        api_bmu_fuel_map=snakemake.params.api_bmu_fuel_map,
    )

    df_bod = get_year_bod(
        base_url,
        filter_bmu_units,
        bmu_carrier_map,
        year=int(snakemake.wildcards.bod_year),
    )

    df_bod.to_csv(snakemake.output.csv)
    logger.info(
        f"Saved Bid/Offer average cost data for the year {snakemake.wildcards.bod_year} to {snakemake.output.csv}"
    )
