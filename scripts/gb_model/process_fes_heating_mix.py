# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Heating technology mix calculator.

This script extracts the heating technology mix from FES workbook and calculates the share of each technology in the final demand.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def process_fes_heatmix(
    fes_data_heatmix_path: str,
    electrified_heating_technologies: dict[str, list[str]],
    scenario: str,
    year_range: list[int],
) -> pd.DataFrame:
    """
    Process heating technology mix

    Args:
        fes_data_heatmix_path(str): Filepath to the heating mix CSV file for each sector
        electrified_heating_technologies (list[str]): List of technologies that contribute to electrified heating demand
        scenario (str): FES scenario considered for the modelling
        year (int): Modeling year
        busmap_path(str): Filepath to the CSV file containing the mapping between the clustered buses and buses in the base PyPSA network

    Returns:
        pd.DataFrame : pandas dataframe containing the share of heating technologies that contribute to electrified heating demand
    """

    # Read the FES data
    fes_data = pd.read_csv(
        fes_data_heatmix_path, index_col=["type", "2020", "scenario_2050"]
    )
    mask = fes_data.index.get_level_values("scenario_2050").str.contains(
        scenario, case=False
    )
    # Filter the data
    fes_scenario_data = (
        fes_data.loc[mask]
        .rename(columns={"data": "2050"})
        .reset_index("2020")
        .reset_index("scenario_2050", drop=True)
        .replace(regex=r"\s*\-\s*", value=float("nan"))
        .astype(float)
    )

    fes_data_grouped = pd.DataFrame(
        float("nan"),
        index=electrified_heating_technologies.keys(),
        columns=fes_scenario_data.columns,
    )
    for tech, to_group in electrified_heating_technologies.items():
        fes_data_grouped.loc[tech] = fes_scenario_data.loc[to_group].sum()

    fes_data_grouped.columns = fes_data_grouped.columns.astype(int).rename("year")
    fes_data_all_years = fes_data_grouped.reindex(
        columns=range(2020, 2051)
    ).interpolate(axis=1)
    fes_data_year_range = fes_data_all_years.loc[:, slice(*year_range)].stack()
    fes_data_share = fes_data_year_range.div(
        fes_data_year_range.groupby("year").sum(), level="year"
    ).rename("share")
    return fes_data_share


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem, year=2022)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    electrified_heating_technologies = snakemake.params.electrified_heating_technologies

    residential_share = process_fes_heatmix(
        fes_data_heatmix_path=snakemake.input.fes_residential_heatmix,
        electrified_heating_technologies=electrified_heating_technologies,
        scenario=snakemake.params.scenario,
        year_range=snakemake.params.year_range,
    )

    services_share = process_fes_heatmix(
        fes_data_heatmix_path=snakemake.input.fes_services_heatmix,
        electrified_heating_technologies=electrified_heating_technologies,
        scenario=snakemake.params.scenario,
        year_range=snakemake.params.year_range,
    )
    heating_mix = pd.concat(
        [residential_share, services_share],
        keys=["residential", "services"],
        names=["sector", "technology", "year"],
    )
    assert np.allclose(heating_mix.groupby(["sector", "year"]).sum(), 1), (
        "Heating mix shares do not sum to 1"
    )
    heating_mix.to_csv(snakemake.output.csv)
