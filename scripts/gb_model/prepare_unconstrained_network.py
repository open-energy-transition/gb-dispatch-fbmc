# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Prepare network for constrained optimization.
"""

import logging
from pathlib import Path

import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

LARGE_NUMBER = 1e6


def _get_intra_gb(df: pd.DataFrame) -> pd.Series:
    return df.bus0.str.startswith("GB") & df.bus1.str.startswith("GB")


def copperplate_gb(network: pypsa.Network) -> None:
    """
    Modify the network in-place to represent a copperplate model of Great Britain by
    setting the capacities of all transmission lines to a very high value.

    Parameters
    ----------
    network : pypsa.Network
        The input PyPSA network representing the GB power system.

    Returns
    -------
    pypsa.Network
        The modified PyPSA network with unconstrained transmission lines.
    """
    network.lines.loc[_get_intra_gb(network.lines), "s_nom"] = LARGE_NUMBER
    network.links.loc[
        _get_intra_gb(network.links) & network.links.carrier.isin(["DC"]), "p_nom"
    ] = LARGE_NUMBER
    logger.info(
        "Set all intra-GB transmission line capacities and DC link capacities to very large values for copperplate model."
    )


def unconstrain_marginal_plant_max_dispatch(network: pypsa.Network) -> None:
    """
    Modify the network in-place to allow the marginal plants in Europe to operate at any capacity,
    so that our "load shedding" equivalent is represented by these plants increasing their output.

    Parameters
    ----------
    network : pypsa.Network
        The input PyPSA network representing the GB power system.

    """
    # Remove GB load shedding plants; we shouldn't need them if we are unconstraining
    # the marginal plants in Europe
    network.remove(
        "Generator",
        network.generators.filter(regex=r"GB \d+ Load Shedding", axis=0).index,
    )
    max_marginal_plant_cost = network.generators[
        (network.generators.set == "PP") & ~network.generators.bus.str.startswith("GB")
    ].marginal_cost.max()
    network.generators.loc[
        network.generators.index.str.contains("Load Shedding"), "marginal_cost"
    ] = max_marginal_plant_cost + 1.0


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem, year=2022)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load network
    network = pypsa.Network(snakemake.input.network)

    copperplate_gb(network)
    if snakemake.params["unconstrain_marginal_eur_plants"]:
        unconstrain_marginal_plant_max_dispatch(network)
    network.export_to_netcdf(snakemake.output.network)
