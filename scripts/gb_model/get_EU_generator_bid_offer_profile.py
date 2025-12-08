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


def filter_interconnectors(df):
    """
    Filter to obtain links between GB and EU
    """
    m1 = df["bus0"].str.startswith("GB")
    m2 = df["bus1"].str.startswith("GB")

    return df[(m1 & ~m2) | (~m1 & m2)].query("carrier == 'DC'")


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
        bus0_marginal_price = marginal_price_profile[connector.bus0]
        bus1_marginal_price = marginal_price_profile[connector.bus1]
        fee = (bus0_marginal_price - bus1_marginal_price).abs()
        fee_profile[idx] = fee

    logger.info("Calculated interconnector fee")
    return fee_profile


def calc_bid_offer_multiplier(
    gb_marginal_electricity_price: float,
    generator_marginal_prices: pd.DataFrame,
    bid_multiplier: dict[str, list],
    offer_multiplier: dict[str, list],
    strike_prices: dict[str, list],
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
    strike_prices: dict[str, list]
        Strike prices for renewable carriers
    """

    # Calculate marginal generator at the EU node
    marginal_gen = (
        (generator_marginal_prices - gb_marginal_electricity_price).abs().sort_values()
    )
    marginal_carrier = marginal_gen.index[0]
    marginal_gen_cost = marginal_gen.iloc[0]

    # Adjust marginal cost for both offer and bid
    if marginal_carrier in strike_prices.keys():
        return strike_prices[marginal_carrier], strike_prices[marginal_carrier]
    elif marginal_carrier in bid_multiplier.keys():
        return marginal_gen_cost * bid_multiplier[
            marginal_carrier
        ], marginal_gen_cost * offer_multiplier[marginal_carrier]
    else:
        return 1, 1


def get_EU_marginal_generator(
    countries: list[str],
    marginal_price_profile: pd.DataFrame,
    unconstrained_result: pypsa.Network,
    bids_and_offers_multipliers: dict[dict[str, list]],
    strike_prices: dict[str, list],
) -> pd.DataFrame:
    countries.remove("GB")
    EU_buses = unconstrained_result.buses.query(
        "country in @countries and carrier == 'AC'"
    ).index
    columns = [f"{x} bid" for x in EU_buses] + [f"{x} offer" for x in EU_buses]

    EU_marginal_gen_profile = pd.DataFrame(
        index=unconstrained_result.snapshots, columns=columns
    )
    bid_multiplier = bids_and_offers_multipliers["bid_multiplier"]
    offer_multiplier = bids_and_offers_multipliers["offer_multiplier"]

    # Calculate bid and offer profile for the marginal generator at each EU node
    for EU_bus in EU_buses:
        generator_marginal_prices = (
            unconstrained_result.generators.query("bus == @EU_bus")
            .groupby("carrier")
            .marginal_cost.mean()
            .round(3)
        )
        EU_marginal_gen_profile[[f"{EU_bus} bid", f"{EU_bus} offer"]] = (
            marginal_price_profile[EU_bus]
            .apply(
                lambda x: calc_bid_offer_multiplier(
                    np.round(x, 3),
                    generator_marginal_prices,
                    bid_multiplier,
                    offer_multiplier,
                    strike_prices,
                )
            )
            .tolist()
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
    strike_prices = pd.read_csv(
        snakemake.input.strike_prices, index_col="carrier"
    ).to_dict()["strike_price_GBP_per_MWh"]
    bids_and_offers_multipliers = snakemake.params.bids_and_offers
    countries = snakemake.params.countries

    marginal_price_profile = extract_marginal_price_profiles(unconstrained_result)

    interconnector_fee_profile = compute_interconnector_fee(
        marginal_price_profile, unconstrained_result
    )

    EU_marginal_gen_profile = get_EU_marginal_generator(
        countries,
        marginal_price_profile,
        unconstrained_result,
        bids_and_offers_multipliers,
        strike_prices,
    )

    interconnector_fee_profile.to_csv(snakemake.output.interconnector_fee)
    logger.info(
        f"Exported interconnector fee profile to {snakemake.output.interconnector_fee}"
    )

    EU_marginal_gen_profile.to_csv(snakemake.output.generator_csv)
    logger.info(
        f"Exported EU marginal generator bid/offer profile to {snakemake.output.generator_csv}"
    )
