# SPDX-FileCopyrightText: gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT

"""
Historical electrified heat demand profile processor.
"""

import logging
from pathlib import Path

import pandas as pd
import xarray as xr

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def get_net_electrified_heat_demand_profile(
    heat_demand: xr.Dataset,
    energy_totals: pd.DataFrame,
    resistive_heater_demand: pd.DataFrame,
) -> pd.DataFrame:
    """
    Get the net electrified heat demand profile by removing future resistive heater demand from historical electrified heat demand.

    If the result is negative, there is more resistive heater demand in future compared to historical.
    When this is subsequently subtracted from the baseline electricity demand, it would then _increase_ the contribution of resistive heating compared to its historical contribution.

    Args:
        heat_demand (xr.Dataset): Profile of hourly heat demand based on heating degree days.
        energy_totals (pd.DataFrame): Historical annual energy demand totals.
        resistive_heater_demand (pd.DataFrame): Future resistive heat demand.

    Returns:
        pd.DataFrame: _description_
    """
    heat_demand_shape = heat_demand / heat_demand.sum("snapshots")

    heat_demand_annual = (
        energy_totals.rename(columns=lambda x: x.removeprefix("electricity").strip())
        .subtract(resistive_heater_demand)
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
    heat_demand = xr.open_dataset(snakemake.input.heat_demand_shape).rename(
        {"node": "bus"}
    )

    energy_totals = (
        pd.read_csv(snakemake.input.energy_totals, index_col="name")
        .mul(1e6)  # TWh to MWh
        .rename_axis(index="bus")
    )
    resistive_heater_demand: dict = {}
    for sector in ["residential", "services"]:
        water_to_space_ratio = energy_totals[f"electricity {sector} water"] / (
            energy_totals[f"electricity {sector} space"]
            + energy_totals[f"electricity {sector} water"]
        )
        df = (
            pd.read_csv(
                snakemake.input[f"{sector}_heat_techs_consumption"],
                index_col=["bus", "year", "technology"],
            )
            .xs(
                (int(snakemake.wildcards.year), "resistive"),
                level=("year", "technology"),
            )
            .squeeze()
        )
        resistive_heater_demand[f"{sector} water"] = df * water_to_space_ratio
        resistive_heater_demand[f"{sector} space"] = df * (1 - water_to_space_ratio)

    net_heat_demand = get_net_electrified_heat_demand_profile(
        heat_demand, energy_totals, pd.DataFrame(resistive_heater_demand)
    )

    net_heat_demand.rename_axis(index="time").to_csv(snakemake.output.csv)
