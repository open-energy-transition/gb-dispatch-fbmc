# SPDX-FileCopyrightText: gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT

"""
EV demand profile processor.

This script prepares regional EV demand profiles.
"""

import logging
from pathlib import Path

import pandas as pd
import xarray as xr

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

HEAT_DEMAND_MAPPING = {"historical": "electricity", "future_resistive": "total"}


def prepare_demand_shape(
    electricity_demand: pd.DataFrame,
    total_heat_demand: pd.DataFrame,
) -> pd.DataFrame:
    """
    Parse and prepare different demand profiles.

    This function processes regional demand data and calculates
    normalized demand profiles (shapes) for use in EV demand modeling.

    Args:
        demand_path (str): Path to the CSV file containing regional
                                demand data with buses as columns
                                   and time periods as index.

    Returns:
        pd.DataFrame: Normalized demand profiles with the same structure
                     as input but with values representing demand shares/proportions.
                     Each column (region) sums to 1.0 across all time periods.

    Processing steps:
        1. Load regional demand data from CSV file
        2. Calculate normalized demand profiles by dividing each region's demand
           by its total annual demand
        3. Return demand shape profiles for use in demand modeling
    """
    net_electricity_demand = electricity_demand.subtract(
        total_heat_demand, fill_value=0
    )
    perc_negative = (
        (net_electricity_demand < 0).sum().sum() / net_electricity_demand.size * 100
    )

    logger.warning(
        f"{perc_negative:.2f}% negative net electricity demand values set to zero"
    )
    demand_shape = (
        net_electricity_demand.clip(lower=0).apply(lambda x: x / x.sum()).fillna(0)
    )

    return demand_shape


def get_net_electrified_heat_demand_profile(
    heat_demand_shape_path: str,
    energy_totals_path: str,
    resistive_heat_demand_path: str,
    year: int,
) -> pd.DataFrame:
    resistive_heat_demand = (
        pd.read_csv(resistive_heat_demand_path, index_col=["sector", "bus", "year"])
        .squeeze()
        .xs(year, level="year")
        .unstack("sector")
        .rename_axis(index="node")
    )
    for use_type in ["space", "water"]:
        resistive_heat_demand = resistive_heat_demand.assign(
            **resistive_heat_demand.rename(columns=lambda x: f"{x} {use_type}")
        )
    heat_demand = xr.open_dataset(heat_demand_shape_path)
    heat_demand_shape = heat_demand / heat_demand.sum("snapshots")
    df_totals = (
        pd.read_csv(energy_totals_path, index_col="name")
        .mul(1e6)
        .rename_axis(index="node")
    )  # TWh to MWh

    heat_demand_annual = (
        df_totals.rename(columns=lambda x: x.removeprefix("electricity").strip())
        .subtract(resistive_heat_demand)
        .to_xarray()
    )
    heat_demand_profile = (
        (heat_demand_shape * heat_demand_annual).to_array().sum("variable").to_pandas()
    )

    return heat_demand_profile


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem, clusters="clustered")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the file path
    historical_demand = pd.read_csv(
        snakemake.input.historical_profile, index_col="time", parse_dates=True
    )

    historical_heat_demand = get_net_electrified_heat_demand_profile(
        snakemake.input.heat_demand_shape,
        snakemake.input.energy_totals,
        snakemake.input.resistive_heat_demand,
        int(snakemake.wildcards.year),
    )

    # Prepare profile shape
    demand_shape = prepare_demand_shape(historical_demand, historical_heat_demand)
    # Save the transport demand profiles
    demand_shape.to_csv(snakemake.output.demand_shape)
    logger.info(
        f"Transport demand profile shapes saved to {snakemake.output.demand_shape}"
    )
