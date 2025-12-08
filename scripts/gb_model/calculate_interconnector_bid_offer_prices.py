# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Calculate interconnector bids and offer prices
"""

import logging
from pathlib import Path

import pandas as pd
import pypsa
import numpy as np

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

def extract_marginal_price_profiles(network):

    ac_buses = network.buses.query("carrier == 'AC'").index
    marginal_price_profile = network.buses_t.marginal_price[ac_buses]

    return marginal_price_profile


def filter_interconnectors(df):
    m1 = df["bus0"].str.startswith("GB")
    m2 = df["bus1"].str.startswith("GB")

    return df[(m1 & ~m2) | (~m1 & m2)].query("carrier == 'DC'")

def compute_interconnector_fee(marginal_price_profile, unconstrained_result):

    interconnectors = filter_interconnectors(unconstrained_result.links)
    fee_profile=pd.DataFrame(index=unconstrained_result.snapshots, columns=interconnectors.index)

    for idx, connector in interconnectors.iterrows():
        bus0_marginal_price = marginal_price_profile[connector.bus0]
        bus1_marginal_price = marginal_price_profile[connector.bus1]
        fee = (bus0_marginal_price - bus1_marginal_price).abs()
        fee_profile[idx] = fee

    return fee_profile


def get_price_setting_carrier(x, generator_marginal_prices, bid_multiplier, offer_multiplier, strike_prices):
    marginal_gen = (generator_marginal_prices - x).abs().sort_values()
    marginal_carrier = marginal_gen.index[0]
    marginal_gen_cost = marginal_gen.iloc[0]
    if marginal_carrier in strike_prices.keys():
        return strike_prices[marginal_carrier], strike_prices[marginal_carrier]
    elif marginal_carrier in bid_multiplier.keys():
        return marginal_gen_cost * bid_multiplier[marginal_carrier], marginal_gen_cost * offer_multiplier[marginal_carrier]
    else:
        return 1,1


def get_EU_marginal_generator(countries, marginal_price_profile, unconstrained_result, bids_and_offers_multipliers, strike_prices):
    countries.remove('GB')
    EU_buses = unconstrained_result.buses.query("country in @countries and carrier == 'AC'").index
    columns = [f"{x} bid" for x in EU_buses] + [f"{x} offer" for x in EU_buses]

    EU_marginal_gen_profile=pd.DataFrame(index=unconstrained_result.snapshots, columns=columns)
    bid_multiplier = bids_and_offers_multipliers['bid_multiplier']
    offer_multiplier = bids_and_offers_multipliers['offer_multiplier']

    for EU_bus in EU_buses:
        generator_marginal_prices=(unconstrained_result.generators
                                .query("bus == @EU_bus")
                                .groupby("carrier")
                                .marginal_cost.mean()
                                .round(3)
                                )   
        EU_marginal_gen_profile[[f"{EU_bus} bid",f"{EU_bus} offer"]]=marginal_price_profile[EU_bus].apply(lambda x: get_price_setting_carrier(np.round(x,3), generator_marginal_prices, bid_multiplier, offer_multiplier, strike_prices)).tolist()
    return EU_marginal_gen_profile

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
    breakpoint()
    df_multiplier = df.map(multiplier).fillna(1)

    df["marginal_cost"] *= df["multiplier"]

    df["marginal_cost"] = (
        df["carrier"]
        .map(strike_prices)
        .where(df["carrier"].isin(strike_prices), df["marginal_cost"])
    )

    return df

def calc_generator_bids_offers(EU_marginal_gen_profile, countries, bids_and_offers_multipliers, strike_prices):
    columns = [f"{x} bid" for x in countries] + [f"{x} offer" for x in countries]
    EU_bid_offer_profile = pd.DataFrame(index=EU_marginal_gen_profile.index, columns=columns)
    for country in countries:
        EU_bid_offer_profile[f"{country} bid"] = _apply_multiplier_and_strike_prices(EU_marginal_gen_profile[country], bids_and_offers_multipliers['bid_multiplier'], strike_prices)



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

    interconnector_fee_profile = compute_interconnector_fee(marginal_price_profile, unconstrained_result)

    EU_marginal_gen_profile = get_EU_marginal_generator(countries, marginal_price_profile, unconstrained_result, bids_and_offers_multipliers, strike_prices)

    EU_marginal_gen_profile.to_csv(snakemake.output.csv)
    logger.info(f"Exported profile to {snakemake.output.csv}")
