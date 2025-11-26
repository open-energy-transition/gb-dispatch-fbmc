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
from scripts.gb_model._helpers import get_lines

logger = logging.getLogger(__name__)


def _prune_lines(
    n: pypsa.Network,
    prune_lines: list[dict[str, int]],
) -> None:
    """
    Prune lines between bus0 and bus1 by setting their s_max_pu to 0.

    Args:
        n (pypsa.Network): The PyPSA network to modify.
        prune_lines (list[dict[str, int]]): The lines to prune.
    """
    # Prune specified lines
    for line in prune_lines:
        mask = get_lines(n.lines, line["bus0"], line["bus1"])
        if mask.any():
            n.lines.loc[mask, "s_max_pu"] = 0
            logger.info(
                f"Pruned line between bus {line['bus0']} and bus {line['bus1']}"
            )
        else:
            logger.warning(
                f"No line found to prune between bus {line['bus0']} and bus {line['bus1']}"
            )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem, year=2022)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load parameters
    prune_lines = snakemake.params.prune_lines

    # Load network
    network = pypsa.Network(snakemake.input.network)

    # Prune lines
    _prune_lines(network, prune_lines)

    network.export_to_netcdf(snakemake.output.network)
