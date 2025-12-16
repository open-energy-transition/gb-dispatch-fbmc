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
from scripts.gb_model._helpers import filter_interconnectors, marginal_costs_bus

logger = logging.getLogger(__name__)


def calculate_interconnector_loss(
    unconstrained_result: pypsa.Network, interconnectors: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculate capacity weighted average of interconnector losses at each foreign market

    Parameters
    ----------
    unconstrained_result: pypsa.Network
        Result of the unconstrained optimization
    interconnectors: pd.DataFrame
        Dataframe of interconnectors and their pypsa parameters
    """
    interconnector_grouping = (
        interconnectors.index.to_series().groupby(interconnectors["bus1"]).agg(list)
    )
    loss_profile = pd.DataFrame(
        index=unconstrained_result.snapshots, columns=interconnector_grouping.index
    )
    for bus in interconnector_grouping.index:
        group = interconnector_grouping.loc[bus]
        # Capacity weighted average of interconnectors connected to each eur node
        loss_profile[bus] = (
            (
                unconstrained_result.links_t.p1[group]
                - unconstrained_result.links_t.p0[group] * -1
            )
            .abs()  # interconnector loss between GB and eur node
            .mul(interconnectors.loc[group].p_nom)
            .sum(axis=1)
            .div(interconnectors.loc[group].p_nom.sum())
        )

    if loss_profile.isna().any().any():
        logger.info(
            f"NaNs found in the loss profile at buses {loss_profile.columns[loss_profile.isna().any()].tolist()}"
        )
        loss_profile = loss_profile.fillna(0)

    logger.info("Calculated the interconnector loss profile")
    return loss_profile


def calc_interconnector_bids_and_offers(
    interconnectors: pd.DataFrame,
    eur_marginal_gen_profile: pd.DataFrame,
    interconnector_fee_profile: pd.DataFrame,
    loss_profile: pd.DataFrame,
    load_shedding_price: float,
) -> pd.DataFrame:
    """
    Calculate the bids and offers for import, export and floating condition of each interconnector

    Parameters
    ----------
        interconnectors: pd.DataFrame
            Dataframe of interconnectors from the pypsa network.links
        eur_marginal_gen_profile: pd.DataFrame
            bid / offer for the marginal generator at each eur bus
        interconnector_fee_profile: pd.DataFrame
            Spread in marginal prices for each interconnector
        loss_profile: pd.DataFrame
            Weighted capacity average of losses of interconnectors connected to each eur bus
        load_shedding_price: float
            Cost of load shedding
    """
    profile_dict = {}

    for connector, row in interconnectors.iterrows():
        gb_bus = row["bus0"]
        eur_bus = row["bus1"]

        eur_highest_marginal_cost = marginal_costs_bus(
            eur_bus, unconstrained_result
        ).sort_values(ascending=False)[0]
        GB_highest_marginal_cost = marginal_costs_bus(
            gb_bus, unconstrained_result
        ).sort_values(ascending=False)[0]

        fee = interconnector_fee_profile[connector]
        # Reset bus shadow price if load generator set the LMP
        p_gb = unconstrained_result.buses_t.marginal_price[gb_bus].apply(
            lambda x: GB_highest_marginal_cost if x >= load_shedding_price else x
        )
        p_eur = unconstrained_result.buses_t.marginal_price[eur_bus].apply(
            lambda x: eur_highest_marginal_cost if x >= load_shedding_price else x
        )
        bid = eur_marginal_gen_profile[f"{eur_bus} bid"]
        offer = eur_marginal_gen_profile[f"{eur_bus} offer"]
        loss = loss_profile[eur_bus]

        # Import bid profile
        profile_dict[(connector, "import_bid")] = (
            p_gb - p_eur * (1 + loss) - bid * (1 + loss)
        )

        # Import offer profile / Import float profile
        profile_dict[(connector, "import_offer")] = profile_dict[
            (connector, "float_import")
        ] = fee + offer * (1 + loss)

        # Export float profile / Export bid profile
        profile_dict[(connector, "float_export")] = profile_dict[
            (connector, "export_bid")
        ] = fee - bid * (1 - loss)

        # Export offer profile
        profile_dict[(connector, "export_offer")] = (
            p_eur * (1 - loss) - p_gb + offer * (1 - loss)
        )

    profile_df = pd.DataFrame(profile_dict)

    return profile_df


def assign_bid_offer(
    unconstrained_result: pypsa.Network,
    interconnectors: pd.DataFrame,
    profiles: pd.DataFrame,
) -> pd.DataFrame:
    """
    Calculate bid/offer for each interconnector based on the status of the interconnector at each time step

    Parameters
    ----------
    unconstrained_result: pypsa.Network
        Result of the unconstrained optimization
    interconnectors: pd.DataFrame
        List of interconnectors between EU and GB
    profiles: pd.DataFrame
        six different bid/offer profiles for each interconnector that can be assigned to it based on it's status
    """

    gb_power = unconstrained_result.links_t.p0[interconnectors.index]
    conditions = [
        gb_power < 0,  # interconnector importing
        gb_power == 0,  # interconnector float
        gb_power > 0,  # interconnector exporting
    ]

    def _filter_profiles(df, key):
        df_filtered = df[df.columns[df.columns.get_level_values(1) == key]]
        df_filtered.columns = df_filtered.columns.droplevel(1)
        return df_filtered

    bid_profiles = [
        _filter_profiles(profiles, "import_bid"),
        _filter_profiles(profiles, "float_export"),
        _filter_profiles(profiles, "export_bid"),
    ]
    offer_profiles = [
        _filter_profiles(profiles, "import_offer"),
        _filter_profiles(profiles, "float_import"),
        _filter_profiles(profiles, "export_offer"),
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
    interconnector_fee_profile = pd.read_csv(
        snakemake.input.interconnector_fee_profile,
        index_col="snapshot",
        parse_dates=True,
    )
    eur_marginal_gen_profile = pd.read_csv(
        snakemake.input.eur_marginal_gen_profile, index_col="snapshot", parse_dates=True
    )
    load_shedding_price = snakemake.params.load_shedding_price

    interconnectors = filter_interconnectors(unconstrained_result.links)

    loss_profile = calculate_interconnector_loss(unconstrained_result, interconnectors)

    profile_df = calc_interconnector_bids_and_offers(
        interconnectors,
        eur_marginal_gen_profile,
        interconnector_fee_profile,
        loss_profile,
        load_shedding_price,
    )

    interconnector_profile = assign_bid_offer(
        unconstrained_result, interconnectors, profile_df
    )

    interconnector_profile.to_csv(snakemake.output.bid_offer_profile)
    logger.info(
        f"Exported interconnector bid/offer profiles to {snakemake.output.bid_offer_profile}"
    )
