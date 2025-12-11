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

from scripts._helpers import configure_logging, set_scenario_config
from scripts.gb_model._helpers import filter_interconnectors, marginal_costs_bus

logger = logging.getLogger(__name__)


def _calculate_interconnector_loss(
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
        # Capacity weighted average of interconnectors connected to each EU node
        loss_profile[bus] = (
            (
                unconstrained_result.links_t.p1[group]
                - unconstrained_result.links_t.p0[group] * -1
            )
            .abs()  # interconnector loss between GB and EU node
            .mul(interconnectors.loc[group].p_nom)
            .sum(axis=1)
            .div(interconnectors.loc[group].p_nom.sum())
        )

    logger.info("Calculated the interconnector loss profile")
    return loss_profile


def calc_interconnector_bids_and_offers(
    interconnectors: pd.DataFrame,
    EU_marginal_gen_profile: pd.DataFrame,
    interconnector_fee_profile: pd.DataFrame,
    loss_profile: pd.DataFrame,
) -> tuple[pd.DataFrame]:
    """
    Calculate the bids and offers for import, export and floating condition of each interconnector

    Parameters
    ----------
        interconnectors: pd.DataFrame
            Dataframe of interconnectors from the pypsa network.links
        EU_marginal_gen_profile: pd.DataFrame
            bid / offer for the marginal generator at each EU bus
        interconnector_fee_profile: pd.DataFrame
            Spread in marginal prices for each interconnector
        loss_profile: pd.DataFrame
            Weighted capacity average of losses of interconnectors connected to each EU bus
    """
    (
        import_bid_profile,
        import_offer_profile,
        float_import_profile,
        export_bid_profile,
        export_offer_profile,
        float_export_profile,
    ) = (
        pd.DataFrame(
            index=interconnector_fee_profile.index, columns=interconnectors.index
        )
        for _ in range(6)
    )

    for connector, row in interconnectors.iterrows():
        if "GB" in row["bus0"]:
            gb_bus = row["bus0"]
            EU_bus = row["bus1"]
        else:
            gb_bus = row["bus1"]
            EU_bus = row["bus0"]

        EU_highest_marginal_cost = marginal_costs_bus(
            EU_bus, unconstrained_result
        ).sort_values(ascending=False)[0]
        GB_highest_marginal_cost = marginal_costs_bus(
            gb_bus, unconstrained_result
        ).sort_values(ascending=False)[0]

        fee = interconnector_fee_profile[connector]
        # Reset bus shadow price if load generator set the LMP
        p_gb = unconstrained_result.buses_t.marginal_price[gb_bus].apply(
            lambda x: GB_highest_marginal_cost if x > 1000 else x
        )
        p_eu = unconstrained_result.buses_t.marginal_price[EU_bus].apply(
            lambda x: EU_highest_marginal_cost if x > 1000 else x
        )
        bid = EU_marginal_gen_profile[f"{EU_bus} bid"]
        offer = EU_marginal_gen_profile[f"{EU_bus} offer"]
        loss = loss_profile[EU_bus]

        # Import bid profile
        import_bid_profile[connector] = p_gb - p_eu * (1 + loss) - bid * (1 + loss)

        # Import offer profile / Import float profile
        import_offer_profile[connector], float_import_profile[connector] = [
            fee + offer * (1 + loss)
        ] * 2

        # Export float profile / Export bid profile
        float_export_profile[connector], export_bid_profile[connector] = [
            fee - bid * (1 - loss)
        ] * 2

        # Export offer profile
        export_offer_profile[connector] = p_eu * (1 - loss) - p_gb + offer * (1 - loss)

    return (
        import_bid_profile,
        import_offer_profile,
        float_import_profile,
        float_export_profile,
        export_bid_profile,
        export_offer_profile,
    )


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
    EU_marginal_gen_profile = pd.read_csv(
        snakemake.input.EU_marginal_gen_profile, index_col="snapshot", parse_dates=True
    )

    interconnectors = filter_interconnectors(unconstrained_result.links)

    loss_profile = _calculate_interconnector_loss(unconstrained_result, interconnectors)

    import_bid, import_offer, float_import, float_export, export_bid, export_offer = (
        calc_interconnector_bids_and_offers(
            interconnectors,
            EU_marginal_gen_profile,
            interconnector_fee_profile,
            loss_profile,
        )
    )

    [
        file.to_csv(name)
        for file, name in zip(
            [
                import_bid,
                import_offer,
                float_import,
                float_export,
                export_bid,
                export_offer,
            ],
            snakemake.output,
        )
    ]
    logger.info(f"Exported interconnector bid/offer profiles to {snakemake.output}")
