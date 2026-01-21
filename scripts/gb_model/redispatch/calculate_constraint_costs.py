# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Calculate constraint costs across all redispatch years
"""

import logging
from pathlib import Path

import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

LOAD_SHEDDING_REGEX = "Load Shedding"


def constraint_cost(networks: list[pypsa.Network], extra_years: int) -> float:
    """
    Calculate constraint costs across all redispatch years

    Parameters
    ----------
    networks: list[pypsa.Network]
        List of solved networks for each redispatch year
    extra_years: int
        Number of extra years included in the redispatch model

    Returns
    -------
    float
        Total constraint costs across all redispatch years
    """
    redispatch_components_regex = ".* ramp (up|down)$"
    load_shedding_regex = ".* Load Shedding$"
    year_constraint_costs = {}
    for network in networks:
        year = network.meta["wildcards"]["year"]
        component_constraint_costs = {}
        for component in ["Generator", "StorageUnit", "Link"]:
            bid_offer_cost = network.get_switchable_as_dense(
                component, "marginal_cost"
            ).filter(regex=redispatch_components_regex)
            generation = (
                network.components[component]
                .dynamic["p0" if component == "Link" else "p"]
                .filter(regex=redispatch_components_regex)
            )
            component_constraint_costs[component] = (
                (bid_offer_cost * generation).sum().sum()
            )

        load_shedding_cost = network.get_switchable_as_dense(
            "Generator", "marginal_cost"
        ).filter(regex=load_shedding_regex)
        load_shedding = network.generators_t.p.filter(regex=load_shedding_regex)
        load_shedding_constraint_cost = (load_shedding_cost * load_shedding).sum().sum()
        component_constraint_costs["Load Shedding"] = load_shedding_constraint_cost

        logger.info(
            f"Constraint costs for redispatch year {year}: "
            + ", ".join(f"{k}: {v:,.0f}" for k, v in component_constraint_costs.items())
            + f", Total: {sum(component_constraint_costs.values()):,.0f}"
        )
        year_constraint_costs[year] = sum(component_constraint_costs.values())
    final_year = max(int(i) for i in year_constraint_costs.keys())
    total_constraint_cost = (
        sum(year_constraint_costs.values())
        + year_constraint_costs[str(final_year)] * extra_years
    )
    logger.info(
        "Total constraint costs across all redispatch years "
        f"(including {extra_years} extra years at year {final_year}): {total_constraint_cost:,.0f}"
    )
    return total_constraint_cost


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem, year=2022)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load input networks and parameters
    networks = [pypsa.Network(f) for f in snakemake.input.networks]
    extra_years = snakemake.params.constraint_cost_extra_years

    total_constraint_cost = constraint_cost(networks, extra_years)

    pd.Series({"constraint_cost": total_constraint_cost}).to_csv(
        snakemake.output.csv, header=False
    )
