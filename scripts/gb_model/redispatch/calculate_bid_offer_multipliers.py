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
from scripts.gb_model.generators.assign_costs import (
    _load_costs,
    _load_fes_carbon_costs,
    _load_fes_power_costs,
)

logger = logging.getLogger(__name__)


def calculate_marginal_costs(
    fes_power_costs_path: str,
    fes_carbon_costs_path: str,
    fes_scenario: str,
    tech_costs_path: str,
    costs_config: dict,
    technology_mapping: dict[str, str],
    start_year: int,
    end_year: int,
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
    df_tech = df_tech.join(fes_power_costs.reset_index().set_index(["carrier", "year"]))
    df_tech = df_tech.join(fes_carbon_costs, on="year")
    df_tech = (
        df_tech.reset_index()
        .set_index("carrier")
        .join(costs[["efficiency", "CO2 intensity"]])
    )

    # Fill missing values for fuel and CO2 intensity
    for tech in df_tech.index.unique().tolist():
        for param in ["CO2 intensity", "fuel"]:
            if (
                tech in costs_config["carrier_gap_filling"][param].keys()
                and df_tech.at[tech, param].isna().any()
            ):
                df_tech.at[tech, param] = costs.at[
                    costs_config["carrier_gap_filling"][param][tech], param
                ]

    # Fill any missing values with default characteristics from the config
    df_tech[costs_config["marginal_cost_columns"]] = df_tech[
        costs_config["marginal_cost_columns"]
    ].apply(lambda x: x.fillna(costs_config["default_characteristics"][x.name]["data"]))

    # Calculate marginal costs
    df_tech["marginal_cost"] = (
        df_tech["VOM"]
        + df_tech["fuel"] / df_tech["efficiency"]
        + df_tech["CO2 intensity"] * df_tech["carbon_cost"] / df_tech["efficiency"]
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


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    start_year = int(Path(snakemake.input.bid_offer_data[0]).stem)
    end_year = int(Path(snakemake.input.bid_offer_data[-1]).stem)

    df_cost = calculate_marginal_costs(
        fes_power_costs_path=snakemake.input.fes_power_costs,
        fes_carbon_costs_path=snakemake.input.fes_carbon_costs,
        fes_scenario=snakemake.params.fes_scenario,
        tech_costs_path=snakemake.input.tech_costs,
        costs_config=snakemake.params.costs_config,
        technology_mapping=snakemake.params.technology_mapping,
        start_year=start_year,
        end_year=end_year,
    )

    df_multipliers = calc_bid_offer_multipliers(
        bid_offer_costs_path=snakemake.input.bid_offer_data, df_cost=df_cost
    )

    df_multipliers.to_csv(snakemake.output.csv)
    logger.info(f"Saved bid and offer multipliers to {snakemake.output.csv}")
