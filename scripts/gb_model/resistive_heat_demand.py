# SPDX-FileCopyrightText: gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT

"""
Resistive heat demand profile processor.
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

HEAT_DEMAND_MAPPING = {"residential_heat": "residential", "iandc_heat": "services"}


def get_resistive_demand(
    eur_demand_annual_path: str, heating_mix_path: str, **heat_demand_paths: str
) -> pd.Series:
    """
    Calculate the future resistive heating electricity demand based on annual heat pump heating demand and the mix of future heating technologies.

    Args:
        eur_demand_annual_path (str): European annual heat demand path.
        heating_mix_path (str): Heating mix path.
        **heat_demand_paths (str): Paths to GB heat demand files.

    Returns:
        pd.Series: Future resistive heating demand.
    """
    eur_demand_annual = pd.read_csv(
        eur_demand_annual_path, index_col=["load_type", "bus", "year"]
    )
    all_demands = {}
    for demand_type, sector in HEAT_DEMAND_MAPPING.items():
        gb_demand_annual = pd.read_csv(
            heat_demand_paths[f"gb_{demand_type}_annual"], index_col=["bus", "year"]
        )
        all_demands[sector] = pd.concat(
            [gb_demand_annual, eur_demand_annual.xs(demand_type, level="load_type")]
        )
    annual_demand = pd.concat(
        all_demands.values(), keys=all_demands.keys(), names=["sector", "bus", "year"]
    )
    heating_mix = pd.read_csv(
        heating_mix_path, index_col=["technology", "sector", "year"]
    ).squeeze()
    resistive_share = heating_mix.loc["resistive"]
    rest_share = (
        heating_mix.drop("resistive", level="technology")
        .groupby(["year", "sector"])
        .sum()
    )
    resistive_demand = (resistive_share / rest_share) * annual_demand.squeeze()

    return resistive_demand


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem, clusters="clustered")
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    resistive_demand = get_resistive_demand(
        snakemake.input.eur_demand_annual,
        snakemake.input.heating_mix,
        gb_residential_heat_annual=snakemake.input.gb_residential_heat_annual,
        gb_iandc_heat_annual=snakemake.input.gb_iandc_heat_annual,
    )
    resistive_demand.rename("p_nom").to_csv(snakemake.output.csv)
