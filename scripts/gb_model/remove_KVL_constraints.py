# SPDX-FileCopyrightText: : gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT
import logging
logger = logging.getLogger(__name__)

def remove_KVL_constraints(n, snapshots, snakemake):
    """
    Add custom extra functionality constraints.
    """
    n.model.remove_constraints("Kirchhoff-Voltage-Law")
    logger.info("Remove KVL constraints from the model")