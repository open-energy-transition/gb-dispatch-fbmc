# SPDX-FileCopyrightText: : gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT
import logging

import pandas as pd
import pypsa
from snakemake.script import Snakemake

logger = logging.getLogger(__name__)


def remove_KVL_constraints(n, snapshots, snakemake):
    """
    Add custom extra functionality constraints to remove KVL constraints.
    """
    n.model.remove_constraints("Kirchhoff-Voltage-Law")
    logger.info("Remove KVL constraints from the model")


def custom_constraints(
    n: pypsa.Network,
    snapshots: pd.Index,
    snakemake: Snakemake,
) -> None:
    """
    Apply custom constraints to the PyPSA network.

    Args:
        n (pypsa.Network): The PyPSA network to modify.
        snapshots (pd.Index): The snapshots of the network.
        snakemake (snakemake.Snakemake): The snakemake object for parameters and config.
    """
    # Apply boundary constraints
    remove_KVL_constraints(n, snapshots, snakemake)
