# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Calculate bid and offer multipliers
"""

import logging
from pathlib import Path

import pandas as pd
import numpy as np

from scripts._helpers import configure_logging, set_scenario_config
from scripts.gb_model.generators.assign_costs import _load_fes_carbon_costs, _load_fes_power_costs, _load_costs

logger = logging.getLogger(__name__)


def calc_bid_offer_multipliers(df_bod: pd.DataFrame, srmc: pd.Series, bmu_carrier_map: dict[str, str]) -> pd.DataFrame:
    """
    Calculate bid and offer multipliers based on Bid / Offer Data and Short run marginal costs

    """

    df_bod['carrier'] = df_bod['nationalGridBmUnit'].map(bmu_carrier_map)
    df_bod = df_bod.dropna(subset=['carrier'])
    df_bod['bid_mul'] = df_bod.apply(lambda x: x['bid'] / srmc.loc[x['carrier']], axis=1)
    df_bod['offer_mul'] = df_bod.apply(lambda x: x['offer'] / srmc.loc[x['carrier']], axis=1)
    df_mul = df_bod.groupby('carrier')[['bid_mul', 'offer_mul']].mean()

    return df_mul

def calculate_marginal_costs(
        fes_power_costs_path: str,
        fes_carbon_costs_path: str, 
        fes_scenario: str, 
        tech_costs_path: str, 
        costs_config: dict,
        technology_mapping: dict[str, str],
    ) -> pd.DataFrame:
    """
    Calculate marginal costs for each technology

    Parameters:
    ----------
    fes_power_costs_path: str
        Filepath to FES power costs CSV
    fes_carbon_costs_path: str
        Filepath to FES carbon costs CSV
    fes_scenario: str
        FES Scenario name
    tech_costs_path: str
        Filepath to technology costs CSV
    costs_config: dict
        Configuration dictionary for costs
    technology_mapping: dict[str, str]
        Mapping of technology names
    """

    # Load FES power and carbon costs
    fes_power_costs = _load_fes_power_costs(fes_power_costs_path, fes_scenario)
    fes_carbon_costs = _load_fes_carbon_costs(fes_carbon_costs_path, fes_scenario)
    costs = _load_costs(tech_costs_path, costs_config)

    # Map FES technology names to PyPSA carriers
    power_tech_map = {v:k.split(" ")[0] for k,v in costs_config['fes_VOM_carrier_mapping'].items() if 'CHP' not in k}
    fes_power_costs['carrier'] = (fes_power_costs.index.get_level_values('technology')
                                  .map(power_tech_map)
                                  .str.upper().map(technology_mapping)
                                )
    
    # Create a multi-index dataframe for technologies and years and merge the FES costs data
    multi_index = pd.MultiIndex.from_product([technology_mapping.values(), np.arange(start_year,end_year+1)], names=['carrier', 'year'])
    df_tech = pd.DataFrame(index=multi_index)
    df_tech = df_tech.join(fes_power_costs.reset_index().set_index(['carrier','year']))
    df_tech = df_tech.join(fes_carbon_costs, on='year')
    df_tech = df_tech.reset_index().set_index('carrier').join(costs[['efficiency','CO2 intensity']])

    # Fill missing values for fuel and CO2 intensity
    for tech in df_tech.index.tolist():
        if tech in costs_config['carrier_gap_filling']['CO2 intensity'].keys() and np.isnan(df_tech.at[tech,'CO2 intensity']):
            df_tech.at[tech, 'CO2 intensity'] = df_tech.at[costs_config['carrier_gap_filling']['CO2 intensity'][tech], 'CO2 intensity']
        if tech in costs_config['carrier_gap_filling']['fuel'] and np.isnan(df_tech.at[tech, 'fuel']):
            df_tech.at[tech, 'fuel'] = costs.at[tech, 'fuel']

    # Fill any missing values with default characteristics from the config
    df_tech[costs_config['marginal_cost_columns']] =df_tech[costs_config['marginal_cost_columns']].apply(lambda x: x.fillna(costs_config['default_characteristics'][x.name]['data']))

    # Calculate marginal costs
    df_tech["marginal_cost"] = (
        df_tech["VOM"]
        + df_tech["fuel"] / df_tech["efficiency"]
        + df_tech["CO2 intensity"] * df_tech["carbon_cost"] / df_tech["efficiency"]
    )

    return df_tech


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    fes_year = snakemake.params.fes_year
    start_year = fes_year 
    end_year = fes_year

    df_tech = calculate_marginal_costs(
        fes_power_costs_path=snakemake.input.fes_power_costs,
        fes_carbon_costs_path=snakemake.input.fes_carbon_costs,
        fes_scenario=snakemake.params.fes_scenario,
        tech_costs_path=snakemake.input.tech_costs,
        costs_config=snakemake.params.costs_config,
        technology_mapping=snakemake.params.technology_mapping,
    )

    breakpoint()

