# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Calculate interconnector bids and offer prices
"""

import logging
from pathlib import Path

import pandas as pd
import pypsa
import numpy as np

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

def extract_marginal_price_profiles(network):

    ac_buses = network.buses.query("carrier == 'AC'").index
    marginal_price_profile = network.buses_t.marginal_price[ac_buses]

    return marginal_price_profile


def filter_interconnectors(df):
    m1 = df["bus0"].str.startswith("GB")
    m2 = df["bus1"].str.startswith("GB")

    return df[(m1 & ~m2) | (~m1 & m2)].query("carrier == 'DC'")

def compute_interconnector_fee(marginal_price_profile, unconstrained_result):

    interconnectors = filter_interconnectors(unconstrained_result.links)
    fee_profile=pd.DataFrame(index=unconstrained_result.snapshots, columns=interconnectors.index)

    for idx, connector in interconnectors.iterrows():
        bus0_marginal_price = marginal_price_profile[connector.bus0]
        bus1_marginal_price = marginal_price_profile[connector.bus1]
        fee = (bus0_marginal_price - bus1_marginal_price).abs()
        fee_profile[idx] = fee

    return fee_profile

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem, year=2022)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load input networks and parameters
    unconstrained_result = pypsa.Network(snakemake.input.unconstrained_result)

    marginal_price_profile = extract_marginal_price_profiles(unconstrained_result)

    interconnector_fee_profile = compute_interconnector_fee(marginal_price_profile, unconstrained_result)

    # network.export_to_netcdf(snakemake.output.network)
    # logger.info(f"Exported network to {snakemake.output.network}")
