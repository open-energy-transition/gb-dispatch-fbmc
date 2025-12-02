# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Regional H2 storage data processor.

This script prepares regional disaggregation of H2 storage data.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.gb_model._helpers import get_regional_distribution

logger = logging.getLogger(__name__)


def prepare_regional_h2_demand_data(
    demand: pd.DataFrame, fixed_supply: pd.DataFrame, regional_reference_data: pd.Series
) -> pd.DataFrame:
    """
    Prepare regional disaggregation of H2 data using reference data distribution patterns.

    Args:
        demand (pd.DataFrame): DataFrame containing annual aggregated H2 demand data indexed by year
        fixed_supply (pd.DataFrame): DataFrame containing annual aggregated fixed H2 supply data indexed by year
        regional_reference_data (pd.Series): Series containing regional data with MultiIndex [bus, year]

    Returns:
        pd.Series:
            Series with MultiIndex [bus, year] containing regionally distributed H2 demand in GWh.
            Each region gets a proportional share of the annual demand based on reference distribution patterns.
    """
    net_demand = demand.sum(axis=1) - fixed_supply.sum(axis=1)

    # Calculate regional distribution of grid-connected electrolysis capacities
    electrolysis_distribution = get_regional_distribution(regional_reference_data.p_nom)

    # Apply regional distribution to electricity demand data
    net_demand_regional = net_demand.mul(electrolysis_distribution)
    assert np.isclose(net_demand_regional.sum(), net_demand.sum()), (
        "Regional H2 demand disaggregation error: totals do not match"
    )

    return net_demand_regional.to_frame("p_set")


def prepare_regional_h2_storage_data(
    storage: pd.DataFrame, regional_reference_data: pd.Series
) -> pd.DataFrame:
    """
    Prepare regional disaggregation of H2 storage data using reference data distribution patterns.

    Args:
        storage (pd.DataFrame): DataFrame containing annual aggregated H2 storage data indexed by year
        regional_reference_data (pd.Series): Series containing regional data with MultiIndex [bus, year]

    Returns:
        pd.Series:
            Series with MultiIndex [bus, year] containing regionally distributed H2 storage capacity in GWh.
            Each region gets a proportional share of the annual storage capacity based on reference distribution patterns.
    """

    # Calculate regional distribution of grid-connected electrolysis capacities
    electrolysis_distribution = get_regional_distribution(regional_reference_data.p_nom)

    # Apply regional distribution to storage data
    regional_storage = storage.squeeze().mul(electrolysis_distribution)
    assert np.isclose(regional_storage.sum(), storage.sum()), (
        "Regional H2 storage disaggregation error: totals do not match"
    )

    return regional_storage.to_frame("e_nom")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem, data_type="unmanaged_charging")
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    demand = pd.read_csv(snakemake.input.demand, index_col="year")
    fixed_supply = pd.read_csv(snakemake.input.fixed_supply, index_col="year")
    storage = pd.read_csv(snakemake.input.storage, index_col="year")
    regional_distribution = pd.read_csv(
        snakemake.input.regional_distribution, index_col=["bus", "year"]
    ).squeeze()
    # Prepare regional H2 data

    regional_demand = prepare_regional_h2_demand_data(
        demand,
        fixed_supply,
        regional_distribution,
    )
    regional_storage = prepare_regional_h2_storage_data(
        storage,
        regional_distribution,
    )
    regional_demand.to_csv(snakemake.output.demand)
    regional_storage.assign(set="Store", carrier="hydrogen").to_csv(
        snakemake.output.storage
    )
