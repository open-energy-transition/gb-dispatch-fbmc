# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Custom constraints provider.

Defines custom constraints for the GB model.
"""

import logging

import pandas as pd
import pypsa
from snakemake.script import Snakemake

from scripts.gb_model._helpers import get_lines

logger = logging.getLogger(__name__)


def set_boundary_constraints(
    n: pypsa.Network,
    snapshots: pd.Index,
    snakemake: Snakemake,
) -> None:
    """
    Limit total line flows across each boundary to match ETYS capacities.

    Args:
        n (pypsa.Network): The PyPSA network to add constraints.
        snapshots (pd.Index): The snapshots of the network.
        snakemake (snakemake.Snakemake): The snakemake object for parameters and config.
    """
    # Load ETYS capacities
    etys_capacities = pd.read_csv(snakemake.input.etys_caps, index_col="boundary_name")
    etys_boundaries = snakemake.params.etys_boundaries_to_lines

    # Define Line-s variable (apparent power flow)
    line_s = n.model["Line-s"]

    for boundary, bus_groups in etys_boundaries.items():
        if boundary not in etys_capacities.index:
            logger.warning(
                f"Boundary '{boundary}' not found in ETYS capacities, skipping."
            )
            continue

        capacity_mw = etys_capacities.loc[boundary, "capability_mw"]

        # Get all lines for this boundary
        boundary_lines_mask = pd.Series(False, index=n.lines.index)

        for buses in bus_groups:
            lines_mask = get_lines(n.lines, buses["bus0"], buses["bus1"])
            if not lines_mask.any():
                logger.warning(
                    f"No lines found for boundary '{boundary}' between "
                    f"buses '{buses['bus0']}' and '{buses['bus1']}'"
                )
            boundary_lines_mask = boundary_lines_mask | lines_mask

        boundary_lines = n.lines[boundary_lines_mask].index

        if boundary_lines.empty:
            logger.warning(
                f"No lines found for boundary '{boundary}', skipping constraint."
            )
            continue

        logger.info(
            f"Boundary {boundary}: {len(boundary_lines)} lines, "
            f"capacity={capacity_mw} MW"
        )

        # Get Line-s for boundary lines
        line_s_boundary = line_s.sel(snapshot=snapshots, Line=boundary_lines)

        # Multiply Line-s by s_max_pu of each line to get actual flow limits
        s_max_pu = n.lines.loc[boundary_lines, "s_max_pu"]
        line_s_boundary_scaled = line_s_boundary * s_max_pu

        # Sum across lines to get total flow at the boundary
        lhs = line_s_boundary_scaled.sum("Line")

        # Add bidirectional constraint: total flow ≤ boundary capability
        n.model.add_constraints(
            lhs <= capacity_mw, name=f"etys_boundary_{boundary}_forward"
        )
        n.model.add_constraints(
            lhs >= -capacity_mw, name=f"etys_boundary_{boundary}_backward"
        )

        logger.info(
            f"Added boundary constraint: -{capacity_mw:.0f} <= '{boundary}' <= {capacity_mw:.0f} MW"
        )


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
    set_boundary_constraints(n, snapshots, snakemake)
