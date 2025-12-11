# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Calculate EU marginal generator bids and offer price profiles
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config
from scripts.gb_model._helpers import filter_interconnectors, marginal_costs_bus

logger = logging.getLogger(__name__)


def _extract_marginal_price_profiles(network: pypsa.Network):
    """
    Extract marginal prices at the buses

    Parameters
    ----------
    network: pypsa.Network
        Unconstrained network optimization result
    """
    ac_buses = network.buses.query("carrier == 'AC'").index
    marginal_price_profile = network.buses_t.marginal_price[ac_buses]

    return marginal_price_profile


def compute_interconnector_fee(
    marginal_price_profile: pd.DataFrame, unconstrained_result: pypsa.Network
) -> pd.DataFrame:
    """
    Compute the spread between the marginal costs of GB and EU node

    Parameters
    ----------
    marginal_price_profile: pd.DataFrame
        Dataframe of marginal costs at each bus
    unconstrained_result: pypsa.Network
        pypsa network with the results of the unconstrained optimization
    """

    interconnectors = filter_interconnectors(unconstrained_result.links)
    fee_profile = pd.DataFrame(
        index=unconstrained_result.snapshots, columns=interconnectors.index
    )

    for idx, connector in interconnectors.iterrows():
        # Generators with highest marginal cost to replace the load shedding generator costs
        EU_highest_marginal_cost = marginal_costs_bus(
            connector.bus1, unconstrained_result
        ).sort_values(ascending=False)[0]
        GB_highest_marginal_cost = marginal_costs_bus(
            connector.bus0, unconstrained_result
        ).sort_values(ascending=False)[0]

        EU_bus_price = marginal_price_profile[connector.bus1].apply(
            lambda x: EU_highest_marginal_cost if x > 1000 else x
        )
        GB_bus_price = marginal_price_profile[connector.bus0].apply(
            lambda x: GB_highest_marginal_cost if x > 1000 else x
        )

        fee_profile[idx] = (GB_bus_price - EU_bus_price).abs()

    logger.info("Calculated interconnector fee")
    return fee_profile


def calc_bid_offer_multiplier(
    gb_marginal_electricity_price: float,
    generator_marginal_prices: pd.DataFrame,
    bid_multiplier: dict[str, list],
    offer_multiplier: dict[str, list],
    renewable_payment_profile: pd.Series,
) -> tuple[float]:
    """
    Calculate bid/offer multiplier profile based on the marginal generator at each EU node for every time step

    Parameters
    ----------
    gb_marginal_electricity_price: float
        Marginal cost of electricity at the GB node for the interconnector
    generator_marginal_prices: pd.DataFrame
        Dataframe of marginal costs for each carrier at the EU node
    bid_multiplier: dict[str, list]
        Multiplier for bids for conventional carriers
    offer_multiplier: dict[str, list]
        Multiplier for offers for conventional carriers
    renewable_payment_profile: pd.Series
        Renewable payment profile for a particular timestamp at the EU node
    """

    # Calculate marginal generator at the EU node
    marginal_gen = (
        (generator_marginal_prices - gb_marginal_electricity_price).abs().sort_values()
    )
    marginal_carrier = marginal_gen.index[0]
    marginal_gen_cost = marginal_gen.iloc[0]

    # Resetting the EU shadow bus price if load shedding had occured (using a random high threshold of 1000)
    # The carrier with the highest marginal cost is used to decide the price rates
    if gb_marginal_electricity_price > 1000:
        marginal_gen_cost = generator_marginal_prices[marginal_carrier]

    # Adjust marginal cost for both offer and bid
    if renewable_payment_profile.index.str.contains(marginal_carrier).any():
        cost = renewable_payment_profile[
            renewable_payment_profile.index.str.contains(marginal_carrier)
        ][0]
        return cost, cost
    elif marginal_carrier in bid_multiplier.keys():
        return marginal_gen_cost * bid_multiplier[
            marginal_carrier
        ], marginal_gen_cost * offer_multiplier[marginal_carrier]
    else:
        return marginal_gen_cost, marginal_gen_cost


def get_EU_marginal_generator(
    marginal_price_profile: pd.DataFrame,
    unconstrained_result: pypsa.Network,
    bids_and_offers_multipliers: dict[dict[str, list]],
    renewable_payment_profile: pd.DataFrame,
) -> pd.DataFrame:
    """
    Obtain the marginal generator at EU node at each timestamp

    Parameters
    ----------
    marginal_price_profile: pd.DataFrame
        Bus shadow price profile indexed by time for each EU bus
    unconstrained_result: pypsa.Network
        Unconstrained optimization model result
    bids_and_offers_multipliers: dict[dict[str, list]]
        Bid and offer multiplier to be applied for conventional generation
    renewable_payment_profile: pd.DataFrame
        Renewable payment profile based on CfD contracts adjusted for electricity market price
    """

    # Filter AC buses at EU nodes
    EU_buses = unconstrained_result.buses.query(
        "country!= 'GB' and carrier == 'AC'"
    ).index

    columns = [f"{x} bid" for x in EU_buses] + [f"{x} offer" for x in EU_buses]
    EU_marginal_gen_profile = pd.DataFrame(
        index=unconstrained_result.snapshots, columns=columns
    )

    # Calculate bid and offer profile for the marginal generator at each EU node
    for EU_bus in EU_buses:
        # Get marginal cost of generators present at the EU bus
        generator_marginal_prices = marginal_costs_bus(EU_bus, unconstrained_result)

        EU_marginal_gen_profile[[f"{EU_bus} bid", f"{EU_bus} offer"]] = (
            marginal_price_profile.apply(
                lambda row: calc_bid_offer_multiplier(
                    np.round(row[EU_bus], 3),
                    generator_marginal_prices,
                    bids_and_offers_multipliers["bid_multiplier"],
                    bids_and_offers_multipliers["offer_multiplier"],
                    renewable_payment_profile.filter(like=EU_bus).loc[row.name],
                ),
                axis=1,
            ).tolist()
        )
    return EU_marginal_gen_profile


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem, year=2022)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load input networks and parameters
    unconstrained_result = pypsa.Network(snakemake.input.unconstrained_result)
    renewable_payment_profile = pd.read_csv(
        snakemake.input.renewable_payment_profile,
        index_col="snapshot",
        parse_dates=True,
    )
    bids_and_offers_multipliers = snakemake.params.bids_and_offers

    marginal_price_profile = _extract_marginal_price_profiles(unconstrained_result)

    interconnector_fee_profile = compute_interconnector_fee(
        marginal_price_profile, unconstrained_result
    )

    EU_marginal_gen_profile = get_EU_marginal_generator(
        marginal_price_profile,
        unconstrained_result,
        bids_and_offers_multipliers,
        renewable_payment_profile,
    )

    interconnector_fee_profile.to_csv(snakemake.output.interconnector_fee)
    logger.info(
        f"Exported interconnector fee profile to {snakemake.output.interconnector_fee}"
    )

    EU_marginal_gen_profile.to_csv(snakemake.output.generator_csv)
    logger.info(
        f"Exported EU marginal generator bid/offer profile to {snakemake.output.generator_csv}"
    )
