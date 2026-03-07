# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Calculate bid and offer multipliers
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.gb_model._helpers import get_scenario_name
from scripts.gb_model.generators.assign_costs import (
    _load_costs,
    _load_fes_carbon_costs,
    _load_fes_power_costs,
    calculate_marginal_costs,
)

logger = logging.getLogger(__name__)


def calculate_costs(
    fes_power_costs_path: str,
    fes_carbon_costs_path: str,
    fes_scenario: str,
    tech_costs_path: str,
    costs_config: dict,
    technology_mapping: dict[str, str],
    start_year: int,
    end_year: int,
    historical_gas_cost: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculate marginal costs for each technology

    Parameters
    ----------
    fes_power_costs_path: str
        Filepath to FES power costs CSV
    fes_carbon_costs_path: str
        Filepath to FES carbon costs CSV
    fes_scenario: str
        FES Scenario name
    tech_costs_path: str
        Filepath to technology costs CSV
    costs_config: dict
        Configuration dictionary for costs
    technology_mapping: dict[str, str]
        Mapping of technology names
    """

    # Load FES power and carbon costs
    fes_power_costs = _load_fes_power_costs(fes_power_costs_path, fes_scenario)
    fes_carbon_costs = _load_fes_carbon_costs(fes_carbon_costs_path, fes_scenario)
    costs = _load_costs(tech_costs_path, costs_config)

    # Map FES technology names to PyPSA carriers
    power_tech_map = {
        v: k.split(" ")[0]
        for k, v in costs_config["fes_VOM_carrier_mapping"].items()
        if "CHP" not in k
    }

    fes_power_costs["carrier"] = (
        fes_power_costs.index.get_level_values("technology")
        .map(power_tech_map)
        .str.upper()
        .map(technology_mapping)
    )

    # Create a multi-index dataframe for technologies and years and merge the FES costs data
    multi_index = pd.MultiIndex.from_product(
        [technology_mapping.values(), np.arange(start_year, end_year + 1)],
        names=["carrier", "year"],
    )
    df_tech = pd.DataFrame(index=multi_index)

    # Merge FES power costs
    df_tech = df_tech.join(fes_power_costs.reset_index().set_index(["carrier", "year"]))

    # Merge technology costs data
    df_tech = (
        df_tech.reset_index()
        .set_index("carrier")
        .join(costs[["efficiency", "CO2 intensity"]])
    ).reset_index()

    gas_fuel_keys = [k for k,v in costs_config['carrier_gap_filling']['fuel'].items() if v == 'gas']
    df_tech.loc[df_tech.carrier.isin(gas_fuel_keys),'fuel'] = df_tech.loc[df_tech.carrier.isin(gas_fuel_keys)]['year'].map(historical_gas_cost.to_dict())

    df_tech = calculate_marginal_costs(
        df_tech,
        costs,
        costs_config,
        fes_carbon_costs,
    )

    # Calculate average marginal cost across the years
    df_cost = df_tech.groupby("carrier")["marginal_cost"].mean()

    logger.info("Calculated average marginal costs for each technology.")
    return df_cost


def calc_bid_offer_multipliers(
    bid_offer_costs_path: list[str],
    df_cost: pd.DataFrame,
) -> pd.DataFrame:
    """
    Calculate bid and offer multipliers as a fraction of the average marginal cost for each technology

    Parameters
    ----------
    bid_offer_costs_path: list[str]
        Filepaths of bid and offer costs data for each technology year-wise
    df_costs: pd.DataFrame
        DataFrame containing the average marginal cost for each technology across years
    """

    df_bid_offer = pd.concat(
        [pd.read_csv(path, index_col="carrier") for path in bid_offer_costs_path]
    )

    df_bid_offer_avg = df_bid_offer.groupby(df_bid_offer.index).mean()

    df_multipliers = df_bid_offer_avg.join(df_cost)

    df_multipliers["bid_multiplier"] = (
        df_multipliers["bid"].abs() / df_multipliers["marginal_cost"]
    )

    df_multipliers["offer_multiplier"] = (
        df_multipliers["offer"] / df_multipliers["marginal_cost"]
    )

    logger.info("Calculated bid and offer multipliers for each technology.")

    return df_multipliers[["bid_multiplier", "offer_multiplier"]]

def get_historical_gas_prices(historical_gas_path, years):
    df = pd.read_excel(historical_gas_path, sheet_name="3.2.1 (Annual)",skiprows=10)
    df = df.query('Year in @y', local_dict={'y':years})
    df.set_index('Year',inplace=True)
    df = df['Major power producers: Natural gas (pence per kWh)\n[Note 4]'] 
    df *= 10 # Convert pence per kWh to GBP per MWh

    return df

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    start_year = int(Path(snakemake.input.bid_offer_data[0]).stem)
    end_year = int(Path(snakemake.input.bid_offer_data[-1]).stem)

    bid_offer_years = [int(Path(x).stem) for x in snakemake.input.bid_offer_data]
    historical_gas_cost = get_historical_gas_prices(
        historical_gas_path = snakemake.input.gas_historical_price,
        years = bid_offer_years

    )

    df_cost = calculate_costs(
        fes_power_costs_path=snakemake.input.fes_power_costs,
        fes_carbon_costs_path=snakemake.input.fes_carbon_costs,
        fes_scenario=get_scenario_name(snakemake),
        tech_costs_path=snakemake.input.tech_costs,
        costs_config=snakemake.params.costs_config,
        technology_mapping=snakemake.params.technology_mapping,
        start_year=start_year,
        end_year=end_year,
        historical_gas_cost=historical_gas_cost
    )

    df_multipliers = calc_bid_offer_multipliers(
        bid_offer_costs_path=snakemake.input.bid_offer_data, df_cost=df_cost
    )

    df_multipliers.to_csv(snakemake.output.csv)
    logger.info(f"Saved bid and offer multipliers to {snakemake.output.csv}")
