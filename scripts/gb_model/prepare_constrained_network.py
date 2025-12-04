# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Prepare network for constrained optimization.
"""

import logging
from pathlib import Path

import pypsa

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def fix_dispatch(n, result):
    """
    Fix dispatch of generators and storage units based on the result of unconstrained optimization

    Parameters
    ----------
    n: pypsa.Network
        Network to finalize
    result: pypsa.Network
        Result of the unconstrained optimization
    """

    for comp in result.components:
        if comp.name not in ["Generator", "StorageUnit"]:
            continue
        p_fix = comp.dynamic.p / comp.static.p_nom

        n.components[comp.name].dynamic.p_max_pu = p_fix
        n.components[comp.name].dynamic.p_min_pu = p_fix


def create_up_down_plants(network, result, bids_and_offers):
    """
    Add generators and storage units components that mimic increase / decrease in dispatch

    Parameters
    ----------
    n: pypsa.Network
        Network to finalize
    result: pypsa.Network
        Result of the unconstrained optimization
    """

    for comp in network.components:
        if comp.name not in ["Generator", "StorageUnit"]:
            continue
        g_up = comp.static.copy()
        g_down = comp.static.copy()

        result_component = result.components[comp.name]
        up_limit = (
            result.get_switchable_as_dense(comp.name, "p_max_pu")
            * result_component.static.p_nom
            - result_component.dynamic.p
        ).clip(0) / result_component.static.p_nom
        down_limit = -result_component.dynamic.p / result_component.static.p_nom

        # Add bid and offer multipliers
        bid_multiplier = bids_and_offers["bid_multiplier"]
        offer_multiplier = bids_and_offers["offer_multiplier"]
        g_up = g_up.assign(multiplier=g_up["carrier"].map(bid_multiplier)).fillna(1)
        g_down = g_down.assign(multiplier=g_up["carrier"].map(offer_multiplier)).fillna(
            1
        )
        g_up["marginal_cost"] *= g_up["multiplier"]
        g_down["marginal_cost"] *= g_down["multiplier"]

        # Add generators that can increase dispatch
        network.add(
            comp.name,
            g_up.index,
            suffix=" ramp up",
            p_max_pu=up_limit.loc[:, g_up.index],
            **g_up.drop("p_max_pu", axis=1),
        )

        # Add generators that can decrease dispatch
        network.add(
            comp.name,
            g_down.index,
            suffix=" ramp down",
            p_min_pu=down_limit.loc[:, g_down.index],
            p_max_pu=0,
            **g_down.drop(["p_max_pu", "p_min_pu"], axis=1),
        )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem, year=2022)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load network
    network = pypsa.Network(snakemake.input.network)
    unconstrained_result = pypsa.Network(snakemake.input.unconstrained_result)
    bids_and_offers = snakemake.params.bids_and_offers

    fix_dispatch(network, unconstrained_result)

    create_up_down_plants(network, unconstrained_result, bids_and_offers)

    network.export_to_netcdf(snakemake.output.network)
