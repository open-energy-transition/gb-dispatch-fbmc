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

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def fetch_api_request_data(url: str, retrieval_message: str) -> pd.DataFrame:
    """
    Fetch Bid Offer Data (BOD) from Elexon API for the specified date range.

    Parameters
    ----------
    url: str
        API request URL
    """

    df = pd.DataFrame()

    max_retries = 3
    for i in range(1,max_retries+1):

        try:
            headers = {'User-Agent': 'Mozilla/5.0'}
            
            response = requests.get(url, headers=headers)
            response.raise_for_status()
            if "Acceptance" in retrieval_message:
                json_data = response.json()['data']
            else:
                json_data = response.json()

            df = pd.DataFrame(json_data) 
            
            logger.info(f"Successfully retrieved {retrieval_message}.")
            break # Exit retry loop if successful

        except requests.exceptions.RequestException as e:

            logger.warning(f" Attempt {i}/{max_retries}: Error in fetching {retrieval_message}: {e}")
            if i == max_retries:
                logger.error(f"Max retries reached. Failed to fetch {retrieval_message}.")
        
    return df


def get_year_bod(base_url: str,start_year: int, end_year: int) -> pd.DataFrame:
    """
    Fetch bid / offer data for an entire year by querying day by day

    Parameters:
    ----------
    base_url: str
        Base URL for the Elexon BMRS API
    start_year: int
        Starting year for which to fetch the data
    end_year: int
        Ending year for which to fetch the data
    """
    dfs = []

    start = datetime(start_year, 1, 1)
    end = datetime(end_year, 12, 31)

    current = start

    while current <= end:
        # Retrieving data month by month
        next_month = current + pd.DateOffset(months=1)

        # API request URL for the month
        url = f"{base_url}/datasets/BOD/stream?from={current}&to={next_month}"
    
        # Fetch data for the current month
        df = fetch_api_request_data(url, retrieval_message=f"Bid/Offer data for the dates {current} to {next_month}")

        dfs.append(df)

        # Reset current date to next month
        current = next_month

    return pd.concat(dfs, ignore_index=True)


# def get_acceptances_data(base_url: str,year: int) -> pd.DataFrame:
#     """
#     Fetch bid / offer data for an entire year by querying day by day

#     Parameters:
#     ----------
#     base_url: str
#         Base URL for the Elexon BMRS API
#     year: int
#         Year for which to fetch the data
#     """
#     dfs = []

#     start = datetime(year, 1, 1)
#     end = datetime(year, 12, 31)

#     current = start

#     while current <= end:
#         for settlement_period in range(1,51):
#             # API request URL for the month
#             url = f"{base_url}/balancing/acceptances/all?settlementDate={current.strftime('%Y-%m-%d')}&settlementPeriod={settlement_period}"

#             # Fetch acceptance data for the current day and settlement period
#             df = fetch_api_request_data(url, retrieval_message=f"Acceptance data for date {current} and settlement period {settlement_period}")

#             dfs.append(df)

#         # Reset current date to next month
#         current = current + pd.DateOffset(days=1)

#     return pd.concat(dfs, ignore_index=True)


# def get_balancing_volumes(base_url: str, year: int) -> pd.DataFrame:
#     """
#     Fetch balancing volumes data from Elexon API

#     Parameters
#     ----------

#     base_url: str
#         Base URL for the Elexon BMRS API
#     """

#     dfs = []

#     start = datetime(year, 1, 1)
#     end = datetime(year, 12, 31)

#     current = start

#     while current <= end:
#         # Retrieving data month by month
#         next_month = current + pd.DateOffset(days=1)

#         # API request URL for the month
#         url = f"{base_url}/balancing/nonbm/volumes?from={current.strftime('%Y-%m-%dT%H:%MZ')}&to={next_month.strftime('%Y-%m-%dT%H:%MZ')}"
        
#         # Fetch data for the current month
#         df = fetch_api_request_data(url, retrieval_message=f"Acceptance data for the dates {current} to {next_month}")

#         dfs.append(df)

#         # Reset current date to next month
#         current = next_month

#     return pd.concat(dfs, ignore_index=True)


def get_market_index(base_url: str, start_year: int, end_year: int, data_provider: str) -> pd.DataFrame:
    """
    Fetch balancing volumes data from Elexon API

    Parameters
    ----------

    base_url: str
        Base URL for the Elexon BMRS API
    """

    dfs = []

    start = datetime(start_year, 1, 1)
    end = datetime(end_year, 12, 31)

    current = start

    while current <= end:
        # Retrieving data year by year
        next_year = current + pd.DateOffset(years=1)

        # API request URL for the year
        url = f"{base_url}/datasets/MID/stream?from={current.strftime('%Y-%m-%dT%H:%MZ')}&to={next_year.strftime('%Y-%m-%dT%H:%MZ')}&dataProviders={data_provider}"
        
        # Fetch data for the current year
        df = fetch_api_request_data(url, retrieval_message=f"Market Index data for the dates {current} to {next_year}")

        dfs.append(df)

        # Reset current date to next year
        current = next_year

    return pd.concat(dfs, ignore_index=True)

def calc_bid_offer_multipliers(df_bod: pd.DataFrame, df_market: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate bid and offer multipliers based on Bid / Offer Data and Market Index Data

    """
    df_units = pd.read_excel("/Users/Sermisha/Downloads/BMUFuelType (2).xlsx")
    bmu_carrier_map = dict(df_units[['NESO BMU ID','REG FUEL TYPE']].values)

    # df_bod = pd.read_csv("BOD_data.csv")
    df_bod['carrier'] = df_bod['nationalGridBmUnit'].map(bmu_carrier_map)

    df_market = pd.read_csv("Market_Index.csv", index_col='startTime')
    df_market.set_index('startTime', inplace=True)
    df_market = df_market[pd.to_datetime(df_market.index).year != 2023].sort_index()
    df_market = df_market.query("price != 0 and volume != 0")

    dict_market = df_market['price'].to_dict()
    dict_volume = df_market['volume'].to_dict()

    df_bod['marketprice'] = df_bod['timeFrom'].map(dict_market)
    df_bod['volume'] = df_bod['timeFrom'].map(dict_volume)

    df_bod['bid_mul'] = df_bod.apply(lambda x: x['bid'] / x['marketprice'] , axis=1)
    df_bod['offer_mul'] = df_bod.apply(lambda x: x['offer'] / x['marketprice'], axis=1)

    df_mul = df_bod.groupby('carrier')[['bid_mul', 'offer_mul']].mean()

    return df_mul

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    base_url = "https://data.elexon.co.uk/bmrs/api/v1"

    # bmu_fuel_type = pd.read_excel(snakemake.input.bmu_fuel_type)
    fes_year = snakemake.params.fes_year
    start_year = fes_year - 4
    end_year = fes_year
    df_bod = get_year_bod(base_url, start_year, end_year)
        # df_market = get_market_index(base_url, year, data_provider=snakemake.params.data_provider)

        # df_mul = calc_bid_offer_multipliers(df_bod, df_market)
    breakpoint()