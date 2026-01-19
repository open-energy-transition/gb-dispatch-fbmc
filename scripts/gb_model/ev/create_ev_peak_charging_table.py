# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
EV unmanaged charging demand peak data processor.

This script processes required EV unmanaged charging demand data from the FES workbook.
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.gb_model._helpers import get_regional_distribution, pre_format

logger = logging.getLogger(__name__)


def parse_ev_unmanaged_charging_data(
    unmanaged_charging_sheet_path: str,
    fes_scenario: str,
    year_range: list,
) -> pd.DataFrame:
    """
    Parse the EV unmanaged charging data from FES workbook to obtain charging demand in the required format.

    Args:
        unmanaged_charging_sheet_path (str): Filepath to the unmanaged charging data CSV file containing
                                           EV charging demand data by scenario and year
        fes_scenario (str): FES scenario name to filter the data for
        year_range (list): Two-element list [start_year, end_year] defining the year range

    Returns:
        pd.DataFrame: DataFrame containing EV unmanaged charging demand indexed by year
                     with 'MW' column representing charging demand in MW for the specified
                     scenario and year range.

    Processing steps:
        1. Load and pre-format unmanaged charging data from CSV file
        2. Filter by scenario and year range, then convert GW to MW
    """

    # Load unmanaged charging data
    df_charging_demand = pd.read_csv(unmanaged_charging_sheet_path)

    # Pre_format the dataframe
    df_charging_demand = pre_format(df_charging_demand)

    # Select unmanaged charging demand data
    unmanaged_charging_demand = df_charging_demand[
        df_charging_demand["Data item"].str.lower()
        == "electric vehicle unmanaged peak demand"
    ]

    # Select scenario
    unmanaged_charging_demand = unmanaged_charging_demand[
        unmanaged_charging_demand["Pathway"].str.lower() == fes_scenario.lower()
    ]

    # Select years in the specified range
    unmanaged_charging_demand = unmanaged_charging_demand[
        unmanaged_charging_demand["year"].between(year_range[0], year_range[-1])
    ]

    # Select only required columns
    unmanaged_charging_demand = unmanaged_charging_demand[["year", "data"]].set_index(
        "year"
    )

    # Convert GW to MW
    unmanaged_charging_demand["data"] = unmanaged_charging_demand["data"] * 1000

    # Rename column as MW
    unmanaged_charging_demand = unmanaged_charging_demand.rename(
        columns={"data": "p_nom"}
    )

    return unmanaged_charging_demand


def prepare_regional_ev_data(
    national_data: pd.DataFrame, regional_data: pd.DataFrame
) -> pd.DataFrame:
    """
    Prepare regional disaggregation of EV data using reference data distribution patterns.

    Args:
        input_path (str): Filepath to the EV data CSV file containing
                           annual aggregated data indexed by year
        reference_data_path (str): Filepath to the reference data CSV file containing
                               regional data with MultiIndex [bus, year]

    Returns:
        pd.Series: Series with MultiIndex [bus, year] containing regionally distributed
                  EV storage capacity in GWh. Each region gets a proportional share of
                  the annual storage capacity based on reference distribution patterns.

    Processing steps:
        1. Load annual storage data and regional reference distribution
        2. Calculate regional distribution patterns from reference data
        3. Apply regional distribution to annual storage capacity
    """
    # Get regional distribution
    regional_dist = get_regional_distribution(regional_data)  # Avoid division by zero

    # Fillna values with 0
    regional_dist = regional_dist.fillna(0)

    # Disaggregate EV data regionally
    regional_ev = regional_dist.squeeze() * national_data.squeeze()

    # Keep original column name
    regional_ev = regional_ev.rename(national_data.columns[0])
    return regional_ev


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    log_suffix = "-" + "_".join(snakemake.wildcards) if snakemake.wildcards else ""
    logger = logging.getLogger(Path(__file__).stem + log_suffix)

    # Load the input paths
    unmanaged_charging_sheet_path = snakemake.input.unmanaged_charging_sheet

    # Load parameters
    fes_scenario = snakemake.params.scenario
    year_range = snakemake.params.year_range

    # Parse EV unmanaged charging demand data
    df_unmanaged_demand = parse_ev_unmanaged_charging_data(
        unmanaged_charging_sheet_path,
        fes_scenario,
        year_range,
    )

    # Write unmanaged charging demand dataframe to csv file
    df_unmanaged_demand.to_csv(snakemake.output.csv)
