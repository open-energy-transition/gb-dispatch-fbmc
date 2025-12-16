# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Calculate eur marginal generator bids and offer price profiles
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config
from scripts.gb_model._helpers import filter_interconnectors, marginal_costs_bus

logger = logging.getLogger(__name__)


def extract_marginal_price_profiles(network: pypsa.Network):
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
    Compute the spread between the marginal costs of GB and eur node

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
        #  Replace the load shedding generator costs with the otherwise most expensive generator
        eur_highest_marginal_cost = marginal_costs_bus(
            connector.bus1, unconstrained_result
        ).sort_values(ascending=False)[0]
        GB_highest_marginal_cost = marginal_costs_bus(
            connector.bus0, unconstrained_result
        ).sort_values(ascending=False)[0]

        eur_bus_price = marginal_price_profile[connector.bus1].apply(
            lambda x: eur_highest_marginal_cost if x > 1000 else x
        )
        GB_bus_price = marginal_price_profile[connector.bus0].apply(
            lambda x: GB_highest_marginal_cost if x > 1000 else x
        )

        fee_profile[idx] = (GB_bus_price - eur_bus_price).abs()

    logger.info("Calculated interconnector fee")
    return fee_profile


def get_gen_marginal_cost(
    generator_marginal_prices: pd.DataFrame,
    gb_marginal_electricity_price: float,
    load_shedding_price: float,
) -> tuple[str, float]:
    """
    To identify marginal generator at the EU bus

    Parameters
    ----------
    generator_marginal_prices: pd.DataFrame
        Dataframe of marginal costs for each carrier at the eur node
    gb_marginal_electricity_price: float
        Marginal cost of electricity at the GB node for the interconnector
    load_shedding_price: float
        Cost of load shedding in GBP/MWh
    """

    # Calculate price difference between generator marginal costs and GB marginal electricity price
    price_spread = (
        (generator_marginal_prices - gb_marginal_electricity_price).abs().sort_values()
    )
    marginal_carrier = price_spread.index[0]
    marginal_gen_cost = gb_marginal_electricity_price

    # Resetting the eur shadow bus price if load shedding had occurred (using a random high threshold of 1000)
    # The carrier with the highest marginal cost is used to decide the price rates
    if gb_marginal_electricity_price >= load_shedding_price:
        marginal_gen_cost = generator_marginal_prices[marginal_carrier]

    return marginal_carrier, marginal_gen_cost


def calc_bid_offer_multiplier(
    gb_marginal_electricity_price: float,
    generator_marginal_prices: pd.DataFrame,
    bid_multiplier: dict[str, list],
    offer_multiplier: dict[str, list],
    renewable_payment_profile: pd.Series,
    load_shedding_price: float,
) -> tuple[float]:
    """
    Calculate bid/offer multiplier profile based on the marginal generator at each eur node for every time step

    Parameters
    ----------
    gb_marginal_electricity_price: float
        Marginal cost of electricity at the GB node for the interconnector
    generator_marginal_prices: pd.DataFrame
        Dataframe of marginal costs for each carrier at the eur node
    bid_multiplier: dict[str, list]
        Multiplier for bids for conventional carriers
    offer_multiplier: dict[str, list]
        Multiplier for offers for conventional carriers
    renewable_payment_profile: pd.Series
        Renewable payment profile for a particular timestamp at the eur node
    load_shedding_price: float
        Cost of load shedding in GBP/MWh
    """

    marginal_carrier, marginal_gen_cost = get_gen_marginal_cost(
        generator_marginal_prices, gb_marginal_electricity_price, load_shedding_price
    )

    # Adjust marginal cost for both offer and bid
    if renewable_payment_profile.index.str.contains(marginal_carrier).any():
        bid_cost = offer_cost = renewable_payment_profile[
            renewable_payment_profile.index.str.contains(marginal_carrier)
        ][0]
    elif marginal_carrier in bid_multiplier.keys():
        bid_cost = marginal_gen_cost * bid_multiplier[marginal_carrier]
        offer_cost = marginal_gen_cost * offer_multiplier[marginal_carrier]
    else:
        bid_cost = offer_cost = marginal_gen_cost
    return bid_cost, offer_cost


def get_eur_marginal_generator(
    marginal_price_profile: pd.DataFrame,
    unconstrained_result: pypsa.Network,
    bids_and_offers_multipliers: dict[dict[str, list]],
    renewable_payment_profile: pd.DataFrame,
    load_shedding_price: float,
) -> pd.DataFrame:
    """
    Obtain the marginal generator at eur node at each timestamp

    Parameters
    ----------
    marginal_price_profile: pd.DataFrame
        Bus shadow price profile indexed by time for each eur bus
    unconstrained_result: pypsa.Network
        Unconstrained optimization model result
    bids_and_offers_multipliers: dict[dict[str, list]]
        Bid and offer multiplier to be applied for conventional generation
    renewable_payment_profile: pd.DataFrame
        Renewable payment profile based on CfD contracts adjusted for electricity market price
    load_shedding_price: float
        Cost of load shedding
    """

    # Filter AC buses at eur nodes
    eur_buses = unconstrained_result.buses.query(
        "country!= 'GB' and carrier == 'AC'"
    ).index

    columns = [f"{x} bid" for x in eur_buses] + [f"{x} offer" for x in eur_buses]
    eur_marginal_gen_profile = pd.DataFrame(
        index=unconstrained_result.snapshots, columns=columns
    )

    # Calculate bid and offer profile for the marginal generator at each eur node
    for eur_bus in eur_buses:
        # Get marginal cost of generators present at the eur bus
        generator_marginal_prices = marginal_costs_bus(eur_bus, unconstrained_result)
        eur_marginal_gen_profile[[f"{eur_bus} bid", f"{eur_bus} offer"]] = (
            marginal_price_profile.apply(
                lambda row: calc_bid_offer_multiplier(
                    np.round(row[eur_bus], 3),
                    generator_marginal_prices,
                    bids_and_offers_multipliers["bid_multiplier"],
                    bids_and_offers_multipliers["offer_multiplier"],
                    renewable_payment_profile.filter(regex=rf"^{eur_bus}").loc[
                        row.name
                    ],
                    load_shedding_price,
                ),
                axis=1,
            ).tolist()
        )
    return eur_marginal_gen_profile


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
    load_shedding_price = (
        snakemake.params.load_shedding_price
    )  # replace with voll from PR #154

    marginal_price_profile = extract_marginal_price_profiles(unconstrained_result)

    interconnector_fee_profile = compute_interconnector_fee(
        marginal_price_profile, unconstrained_result
    )

    eur_marginal_gen_profile = get_eur_marginal_generator(
        marginal_price_profile,
        unconstrained_result,
        bids_and_offers_multipliers,
        renewable_payment_profile,
        load_shedding_price,
    )

    interconnector_fee_profile.to_csv(snakemake.output.interconnector_fee)
    logger.info(
        f"Exported interconnector fee profile to {snakemake.output.interconnector_fee}"
    )

    eur_marginal_gen_profile.to_csv(snakemake.output.generator_csv)
    logger.info(
        f"Exported eur marginal generator bid/offer profile to {snakemake.output.generator_csv}"
    )
