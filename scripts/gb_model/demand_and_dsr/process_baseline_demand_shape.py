# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Baseline electricity demand profile processor.
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

HEAT_DEMAND_MAPPING = {"historical": "electricity", "future_resistive": "total"}


def prepare_demand_shape(
    electricity_demand: pd.DataFrame,
    resistive_heater_demand_change: pd.DataFrame,
) -> pd.DataFrame:
    """
    Parse and prepare different demand profiles.

    This function processes regional demand data and calculates the normalized demand profile / shape for baseline electricity demand.

    Args:
        electricity_demand (pd.DataFrame):
            historical electricity demand profile.
        resistive_heater_demand_change (pd.DataFrame):
            Estimated historical resistive heater demand minus future resistive heater demand.
            This is subtracted from historical electricity demand.
            If negative, then it implies more resistive heater demand in future than in the historical year.

    Returns:
        pd.DataFrame:
            Normalized demand profiles with the same structure as input but with values representing demand shares/proportions.
            Each column (region) sums to 1.0 across all time periods.

    Processing steps:
        1. Load regional demand data from CSV file
        2. Calculate normalized demand profiles by dividing each region's demand
           by its total annual demand
        3. Return demand shape profiles for use in demand modeling
    """
    net_electricity_demand = electricity_demand.subtract(
        resistive_heater_demand_change, fill_value=0
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
    resistive_heater_demand_change = pd.read_csv(
        snakemake.input.resistive_heater_demand, index_col="time", parse_dates=True
    )

    demand_shape = prepare_demand_shape(
        historical_demand, resistive_heater_demand_change
    )

    demand_shape.to_csv(snakemake.output.demand_shape)
