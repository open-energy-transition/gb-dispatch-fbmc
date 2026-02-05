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
from scripts.gb_model.generators.assign_costs import _load_fes_carbon_costs, _load_fes_power_costs, _load_costs

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


def get_year_bod(base_url: str, filter_bmu_units:str, start_year: int, end_year: int) -> pd.DataFrame:
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
        url = f"{base_url}/datasets/BOD/stream?from={current}&to={next_month}{filter_bmu_units}"

        # Fetch data for the current month
        df = fetch_api_request_data(url, retrieval_message=f"Bid/Offer data for the dates {current} to {next_month}")

        dfs.append(df)

        # Reset current date to next month
        current = next_month

    return pd.concat(dfs, ignore_index=True)


def calc_bid_offer_multipliers(df_bod: pd.DataFrame, srmc: pd.Series, bmu_carrier_map: dict[str, str]) -> pd.DataFrame:
    """
    Calculate bid and offer multipliers based on Bid / Offer Data and Short run marginal costs

    """

    df_bod['carrier'] = df_bod['nationalGridBmUnit'].map(bmu_carrier_map)
    df_bod = df_bod.dropna(subset=['carrier'])
    df_bod['bid_mul'] = df_bod.apply(lambda x: x['bid'] / srmc.loc[x['carrier']], axis=1)
    df_bod['offer_mul'] = df_bod.apply(lambda x: x['offer'] / srmc.loc[x['carrier']], axis=1)
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
    start_year = fes_year 
    end_year = fes_year

    pypsa_network = pypsa.Network(snakemake.input.compose_network)
    srmc = pd.concat([pypsa_network.generators.groupby('carrier').marginal_cost.mean(), pypsa_network.storage_units.groupby('carrier').marginal_cost.mean()])

    technology_mapping = snakemake.params.technology_mapping

    url = f"{base_url}/reference/bmunits/all"
    df_bmu = fetch_api_request_data(url, retrieval_message="BMU Unit Data")
    df_bmu['carrier'] = df_bmu['fuelType'].map(technology_mapping)
    df_bmu_filtered = df_bmu.loc[df_bmu['carrier'].isin(technology_mapping.values())]
    bmu_units = df_bmu_filtered['nationalGridBmUnit'].unique().tolist()
    filter_bmu_units = ",".join([(f"&bmUnit={x}") for x in bmu_units]).replace(",","")
    bmu_carrier_map = dict(df_bmu_filtered[['nationalGridBmUnit','carrier']].values)

    # df_bod = get_year_bod(base_url, filter_bmu_units, start_year, end_year)
        # df_market = get_market_index(base_url, year, data_provider=snakemake.params.data_provider)
    df_bod = pd.read_csv("BOD_new.csv")
    df_mul = calc_bid_offer_multipliers(df_bod, srmc, bmu_carrier_map)

    fes_power_costs = _load_fes_power_costs(snakemake.input.fes_power_costs, snakemake.params.fes_scenario)
    fes_carbon_costs = _load_fes_carbon_costs(snakemake.input.fes_carbon_costs, snakemake.params.fes_scenario)
    costs = _load_costs(snakemake.input.tech_costs, snakemake.params.costs_config)

    costs_config = snakemake.params.costs_config
    power_tech_map = {v:k.split(" ")[0] for k,v in costs_config['fes_VOM_carrier_mapping'].items() if 'CHP' not in k}
    fes_power_costs['carrier'] = fes_power_costs.index.get_level_values('technology').map(power_tech_map).str.upper().map(technology_mapping)

    multi_index = pd.MultiIndex.from_product([technology_mapping.values(), np.arange(start_year,end_year+1)], names=['carrier', 'year'])
    df_tech = pd.DataFrame(index=multi_index)
    df_tech = df_tech.join(fes_power_costs.reset_index().set_index(['carrier','year']))
    df_tech = df_tech.join(fes_carbon_costs, on='year')
    df_tech = df_tech.reset_index().set_index('carrier').join(costs[['efficiency','CO2 intensity']])

    for tech in df_tech.index.tolist():
        if tech in costs_config['carrier_gap_filling']['CO2 intensity'].keys() and np.isnan(df_tech.at[tech,'CO2 intensity']):
            df_tech.at[tech, 'CO2 intensity'] = df_tech.at[costs_config['carrier_gap_filling']['CO2 intensity'][tech], 'CO2 intensity']
        if tech in costs_config['carrier_gap_filling']['fuel'] and np.isnan(df_tech.at[tech, 'fuel']):
            df_tech.at[tech, 'fuel'] = costs.at[tech, 'fuel']

    df_tech[costs_config['marginal_cost_columns']] =df_tech[costs_config['marginal_cost_columns']].apply(lambda x: x.fillna(costs_config['default_characteristics'][x.name]['data']))
    df_tech["marginal_cost"] = (
        df_tech["VOM"]
        + df_tech["fuel"] / df_tech["efficiency"]
        + df_tech["CO2 intensity"] * df_tech["carbon_cost"] / df_tech["efficiency"]
    )