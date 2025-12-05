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


def fix_dispatch(constrianed_network, unconstrained_result):
    """
    Fix dispatch of generators and storage units based on the result of unconstrained optimization

    Parameters
    ----------
    constrained_network: pypsa.Network
        Constrained network to finalize
    unconstrained_result: pypsa.Network
        Result of the unconstrained optimization
    """

    for comp in unconstrained_result.components:
        if comp.name not in ["Generator", "StorageUnit"]:
            continue
        p_fix = comp.dynamic.p / comp.static.p_nom

        constrianed_network.components[comp.name].dynamic.p_max_pu = p_fix
        constrianed_network.components[comp.name].dynamic.p_min_pu = p_fix
    
    logger.info("Fixed the dispatch of generators")


def _apply_multiplier_and_strike_prices(
    df: pd.DataFrame,
    multiplier: dict[str, float],
    strike_prices: dict[str, float],
) -> pd.DataFrame:
    """
    Apply bid/offer multiplier and strike prices

    Parameters
    ----------
    df: pd.DataFrame
        Generator dataframe
    multiplier: dict[str, float]
        Mapping of conventional carrier to multiplier
    strike_prices: dict[str, float]
        Mapping of renewable carrier to strike price
    """
    df = df.assign(multiplier=df["carrier"].map(multiplier)).fillna(1)

    df["marginal_cost"] *= df["multiplier"]

    df["marginal_cost"] = (
        df["carrier"]
        .map(strike_prices)
        .where(df["carrier"].isin(strike_prices), df["marginal_cost"])
    )

    return df


def create_up_down_plants(
    constrained_network: pypsa.Network,
    unconstrained_result: pypsa.Network,
    bids_and_offers: dict[str, float],
    strike_prices: dict[str, float],
):
    """
    Add generators and storage units components that mimic increase / decrease in dispatch

    Parameters
    ----------
    constrained_network: pypsa.Network
        Constrained network to finalize
    unconstrained_result: pypsa.Network
        Result of the unconstrained optimization
    bids_and_offers: dict[str, float]
        Bid and offer multipliers for conventional carriers
    strike_prices: dict[str, float]
        Strike prices for renewable carriers
    """

    for comp in constrained_network.components:
        if comp.name not in ["Generator", "StorageUnit"]:
            continue
        g_up = comp.static.copy()
        g_down = comp.static.copy()

        # Compute dispatch limits for the up and down generators
        result_component = unconstrained_result.components[comp.name]
        up_limit = (
            unconstrained_result.get_switchable_as_dense(comp.name, "p_max_pu")
            * result_component.static.p_nom
            - result_component.dynamic.p
        ).clip(0) / result_component.static.p_nom
        down_limit = -result_component.dynamic.p / result_component.static.p_nom

        # Add bid/offer multipliers and set strike prices
        g_up = _apply_multiplier_and_strike_prices(
            g_up, bids_and_offers["bid_multiplier"], strike_prices
        )
        g_down = _apply_multiplier_and_strike_prices(
            g_down, bids_and_offers["offer_multiplier"], strike_prices
        )

        # Add generators that can increase dispatch
        constrained_network.add(
            comp.name,
            g_up.index,
            suffix=" ramp up",
            p_max_pu=up_limit.loc[:, g_up.index],
            **g_up.drop("p_max_pu", axis=1),
        )

        # Add generators that can decrease dispatch
        constrained_network.add(
            comp.name,
            g_down.index,
            suffix=" ramp down",
            p_min_pu=down_limit.loc[:, g_down.index],
            p_max_pu=0,
            **g_down.drop(["p_max_pu", "p_min_pu"], axis=1),
        )

        logger.info("Added generators that can mimic increase and decrease in dispatch")

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem, year=2022)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load input networks and parameters
    network = pypsa.Network(snakemake.input.network)
    unconstrained_result = pypsa.Network(snakemake.input.unconstrained_result)
    bids_and_offers = snakemake.params.bids_and_offers
    strike_prices = pd.read_csv(
        snakemake.input.strike_prices, index_col="carrier"
    ).to_dict()["strike_price_GBP_per_MWh"]

    fix_dispatch(network, unconstrained_result)

    create_up_down_plants(network, unconstrained_result, bids_and_offers, strike_prices)

    network.export_to_netcdf(snakemake.output.network)
    logger.info(f"Exported network to {snakemake.output.network}")
