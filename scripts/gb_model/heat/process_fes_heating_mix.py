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
    uptake_curve: pd.DataFrame,
) -> pd.DataFrame:
    """
    Process heating technology mix

    Args:
        fes_data_heatmix_path(str): Filepath to the heating mix CSV file for each sector
        electrified_heating_technologies (list[str]): List of technologies that contribute to electrified heating demand
        scenario (str): FES scenario considered for the modelling
        year (int): Modeling year
        uptake_curve (pd.DataFrame): Residential heat pump uptake curve from FES data, to interpolate heating technology mix over time

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
        index=pd.Index(electrified_heating_technologies.keys(), name="technology"),
        columns=fes_scenario_data.columns,
    )
    for tech, to_group in electrified_heating_technologies.items():
        fes_data_grouped.loc[tech] = fes_scenario_data.loc[to_group].sum()

    fes_data_interpolated = (
        (uptake_curve.to_xarray() * fes_data_grouped.to_xarray())
        .to_array()
        .sum("variable")
    )
    fes_data_year_range = fes_data_interpolated.sel(year=slice(*year_range)).to_series()

    fes_data_share = fes_data_year_range.div(
        fes_data_year_range.groupby("year").sum(), level="year"
    ).rename("share")

    return fes_data_share


def get_uptake_curve(residential_uptake_trend_path: str, scenario: str) -> pd.DataFrame:
    """
    Get the residential heat pump uptake curve from FES data.

    This will be applied to residential and service sectors to interpolate in heating technology mix over time.

    Args:
        residential_uptake_trend_path(str): Filepath to the residential heat pump uptake trend CSV file
        scenario (str): FES scenario considered for the modelling
    Returns:
        pd.Series : pandas series containing the residential heat pump uptake curve indexed by year
    """
    uptake_data = pd.read_csv(residential_uptake_trend_path)
    uptake_curve = (
        uptake_data[uptake_data.scenario.str.lower() == scenario].set_index("year").data
    )
    uptake_curve_normalised = (uptake_curve - uptake_curve.loc[2020]) / (
        uptake_curve.loc[2050] - uptake_curve.loc[2020]
    )
    # we will apply the normalised curve to 2050 data and the 1-uptake curve to 2020 data to get the final quantity of each technology
    dual_curves = pd.concat(
        [uptake_curve_normalised, 1 - uptake_curve_normalised],
        keys=["2050", "2020"],
        axis=1,
    )
    return dual_curves


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    electrified_heating_technologies = snakemake.params.electrified_heating_technologies
    uptake_curve = get_uptake_curve(
        snakemake.input.fes_hp_uptake_trend, snakemake.params.scenario
    )
    residential_share = process_fes_heatmix(
        fes_data_heatmix_path=snakemake.input.fes_residential_heatmix,
        electrified_heating_technologies=electrified_heating_technologies,
        scenario=snakemake.params.scenario,
        year_range=snakemake.params.year_range,
        uptake_curve=uptake_curve,
    )

    services_share = process_fes_heatmix(
        fes_data_heatmix_path=snakemake.input.fes_services_heatmix,
        electrified_heating_technologies=electrified_heating_technologies,
        scenario=snakemake.params.scenario,
        year_range=snakemake.params.year_range,
        uptake_curve=uptake_curve,
    )
    heating_mix = pd.concat(
        [residential_share, services_share],
        keys=["residential", "services"],
        names=["sector", "year", "technology"],
    )
    assert np.allclose(heating_mix.groupby(["sector", "year"]).sum(), 1), (
        "Heating mix shares do not sum to 1"
    )
    heating_mix.to_csv(snakemake.output.csv)
