# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Calculate interconnector bids and offer prices
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config
from scripts.gb_model._helpers import filter_interconnectors

logger = logging.getLogger(__name__)


def assign_bid_offer(
    unconstrained_result: pypsa.Network,
    interconnectors: pd.DataFrame,
    profiles: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """
    Calculate bid/offer for each interconnector based on the status of the interconnector at each time step

    Parameters
    ----------
    unconstrained_result: pypsa.Network
        Result of the unconstrained optimization
    interconnectors: pd.DataFrame
        List of interconnectors between EU and GB
    profiles: dict[str, pd.DataFrame]
        six different bid/offer profiles that can be assigned to the interconnector based on it's status
    """

    gb_power = unconstrained_result.links_t.p0[interconnectors.index]
    conditions = [
        gb_power < 0,  # interconnector importing
        gb_power == 0,  # interconnector float
        gb_power > 0,  # interconnector exporting
    ]
    bid_profiles = [
        profiles["import_bid"],
        profiles["float_export"],
        profiles["export_bid"],
    ]
    offer_profiles = [
        profiles["import_offer"],
        profiles["float_import"],
        profiles["export_offer"],
    ]
    # for connector in interconnectors:
    bid_profiles = pd.DataFrame(
        np.select(conditions, bid_profiles), columns=interconnectors.index
    ).add_suffix(" bid")
    offer_profiles = pd.DataFrame(
        np.select(conditions, offer_profiles), columns=interconnectors.index
    ).add_suffix(" offer")

    interconnector_profile = pd.concat([bid_profiles, offer_profiles], axis=1)
    interconnector_profile.index = unconstrained_result.snapshots
    logger.info(
        "Assigned the bid/offer to each interconnector based on the status of the interconnector"
    )

    return interconnector_profile


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem, year=2022)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load input networks and parameters
    unconstrained_result = pypsa.Network(snakemake.input.unconstrained_result)

    interconnectors = filter_interconnectors(unconstrained_result.links)

    profile_names = [
        "import_bid",
        "import_offer",
        "float_import",
        "float_export",
        "export_bid",
        "export_offer",
    ]

    profiles = {
        f: pd.read_csv(x, index_col="snapshot", parse_dates=True)
        for f, x in zip(profile_names, snakemake.input.bid_offer_profiles)
    }

    interconnector_profile = assign_bid_offer(
        unconstrained_result, interconnectors, profiles
    )

    interconnector_profile.to_csv(snakemake.output.csv)
    logger.info(
        f"Exported CSV of interconnector bid/offer profiles to {snakemake.output.csv}"
    )
