# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT


"""
GSP-level data table generator.

This is a script to combine the BB1 sheet with the BB2 (metadata) sheet of the FES workbook.
"""

import logging
from pathlib import Path

import geopandas as gpd
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.gb_model._helpers import (
    get_scenario_name,
    map_points_to_regions,
    strip_str,
)

logger = logging.getLogger(__name__)


def parse_inputs(
    bb1_path: str,
    bb2_path: str,
    gsp_coordinates_path: str,
    extra_gsp_coordinates: dict,
    manual_gsp_mapping: dict,
    fes_scenario: str,
    year_range: list,
) -> pd.DataFrame:
    """
    Parse the input data to the required format.

    Args:
        bb1_path (str): path of extracted sheet BB1 of the FES workbook
        bb2_path (str): path of extracted sheet BB2 of the FES workbook
        gsp_coordinates_path (str): path of the GSP supply point coordinates file
        fes_scenario (str): FES scenario
    """

    df_bb2 = pd.read_csv(bb2_path)

    # First step: extract the ID numbers from the Parameter column and set it as the index (it is the only unique identifier for table BB2)
    df_bb2 = (
        df_bb2.set_index(
            ["Template", "Technology", "Technology Detail", "Parameter"], append=True
        )
        .squeeze()
        .unstack("Parameter")
    )
    df_bb2_pivoted = (
        df_bb2.bfill()
        .where(~df_bb2["Building Block ID Number"].isnull())
        .dropna(how="all")
        .reset_index()
        .set_index("Building Block ID Number")
        .drop("level_0", axis=1)
        .apply(strip_str)
    )

    df_bb1 = pd.read_csv(bb1_path)
    df_bb1 = df_bb1.apply(strip_str)
    df_bb1_scenario = df_bb1[
        (df_bb1["FES Pathway"].str.lower() == fes_scenario.lower())
        & (df_bb1["year"].between(year_range[0], year_range[1], inclusive="both"))
    ]
    non_data_cols = df_bb1_scenario.columns.drop("data")
    if (duplicates := df_bb1_scenario[non_data_cols].duplicated()).any():
        # Manual inspection suggests these are true duplicates that should be summed
        logger.warning(
            f"There are {duplicates.sum()} duplicate rows in BB1. These will be summed."
        )
    df_bb1_scenario_no_dups = df_bb1_scenario.groupby(
        non_data_cols.tolist(), as_index=False
    )["data"].sum()

    df_bb1_bb2_scenario = pd.merge(
        df_bb1_scenario_no_dups,
        df_bb2_pivoted,
        left_on="Building Block ID Number",
        right_index=True,
    )
    assert len(df_bb1_bb2_scenario) == len(df_bb1_scenario_no_dups), (
        "Some Building Blocks in BB1 are not present in BB2"
    )

    # We allow cases where there is only a partial match ("Number" vs "Number of" by comparing string starts)
    units_match = df_bb1_bb2_scenario.apply(
        lambda x: x.Units.startswith(x.Unit), axis=1
    )
    assert (units_match).all(), (
        "Mapping of building blocks between BB1 and BB2 may be incorrect as some units do not match:\n"
        f"{df_bb1_bb2_scenario[~units_match][['Unit', 'Units']]}"
    )

    df_bb1_bb2_scenario = df_bb1_bb2_scenario.drop(
        columns=["Units", "Building Block ID Number"]
    )

    df_gsp_coordinates = pd.read_csv(gsp_coordinates_path)
    df_gsp_coordinates = df_gsp_coordinates.apply(strip_str)
    extra_gsp_coordinates_df = (
        pd.DataFrame.from_dict(extra_gsp_coordinates, orient="index")
        .rename_axis(index="Name")
        .reset_index()
    )
    df_gsp_coordinates = (
        # There are cases of duplicate GSPs where the lat and lon information is the same but the GSP ID and GSP group are slightly different
        pd.concat([df_gsp_coordinates, extra_gsp_coordinates_df])
        .drop_duplicates(subset=["Name", "Latitude", "Longitude"])
        .dropna(subset=["Latitude", "Longitude"])
        .set_index("Name")
        .reset_index()
    )
    if (dups := df_gsp_coordinates.Name.duplicated()).any():
        logger.error(
            f"There are duplicate GSP names with different lat/lons in the GSP coordinates file:\n{df_gsp_coordinates[dups]}"
        )

    df_bb1_bb2_scenario["GSP"] = df_bb1_bb2_scenario["GSP"].replace(manual_gsp_mapping)

    df_bb1_bb2_with_lat_lon = pd.merge(
        df_bb1_bb2_scenario,
        df_gsp_coordinates,
        left_on="GSP",
        right_on="Name",
    )

    # Missing data checks.
    # We won't raise errors here as we are willing to accept some missing data for now
    missing_lat_lon = df_bb1_bb2_with_lat_lon[
        df_bb1_bb2_with_lat_lon[["Latitude", "Longitude"]].isnull().any(axis=1)
    ].GSP.unique()
    if len(missing_lat_lon) > 0:
        raise ValueError(
            f"The following GSPs are missing latitude and/or longitude information: {missing_lat_lon}.\n"
            "Please update the GSP coordinates file or provide extra coordinates via the `grid-supply_points.fill_lat_lons` configuration option."
        )

    missing_gsps = set(df_bb1_bb2_scenario.GSP).difference(df_bb1_bb2_with_lat_lon.GSP)
    if missing_gsps:
        logger.warning(
            f"The following GSPs are missing from the GSP coordinates file: {missing_gsps}."
            "Their data will be distributed later across other GSPs in the same TO region or across the whole country."
        )
    df_final = pd.concat(
        [df_bb1_bb2_scenario.query("GSP in @missing_gsps"), df_bb1_bb2_with_lat_lon]
    )

    return df_final, df_bb1_scenario_no_dups


def _merge_gsps(df: pd.DataFrame, all_occurences: pd.DataFrame, gsps: str):
    """
    Merge multiple GSPs into a single GSP, by dissolving the geometries

    Parameters
    ----------
    df: pd.DataFrame
        The dataframe containing the GSPs to merge
    all_occurences: pd.DataFrame
        All occurences of the GSPs to merge (i.e. all rows where the GSPs to merge are mentioned in the "GSPs" column)
    gsps: str
        The GSPs to merge, as a single string with "|" separating the different GSPs (e.g. "GSP1|GSP2|GSP3")s
    """

    # Dissolve the rows of GSPs into a single row 
    df.loc[df["GSPs"] == gsps, "geometrycoord"] = (
        df.loc[all_occurences.index]
        .set_geometry("geometrycoord")
        .dissolve()
        .geometrycoord.values[0]
    )
    
    # Assign the GSP name as the concantenation of the GSP names of the merged GSPs, separated by "|"
    df.loc[df["GSPs"] == gsps, "GSP"] = all_occurences[
        all_occurences.GSPs != gsps
    ].GSP.str.cat(sep="|")

    # Drop the other rows of the merged GSPs, keeping only the dissolved row
    indices = all_occurences.index.tolist()
    indices = [
        x
        for x in all_occurences.index
        if x not in df.loc[df["GSPs"] == gsps].index.tolist()
    ]

    df.drop(index=indices, inplace=True)

    return df


def create_gsp_shapefile(
    gsp_coordinates_path: str,
    gsp_shapes_path: str,
    df_bb1: pd.DataFrame,
    regions: gpd.GeoDataFrame,
    gsp_mapping: dict,
    extra_gsp_coordinates: dict,
    combine_busbars: dict,
):
    df_gsp_coordinates = pd.read_csv(gsp_coordinates_path)
    gdf_gsps = gpd.GeoDataFrame(
        df_gsp_coordinates,
        geometry=gpd.points_from_xy(
            df_gsp_coordinates.Longitude, df_gsp_coordinates.Latitude
        ),
        crs="EPSG:4326",
    )

    df_gsp_shapes = gpd.read_file(gsp_shapes_path)

    df_bb1_gsp = pd.DataFrame(data=df_bb1.GSP.unique(), columns=["GSP"])
    df_bb1_gsp["GSP"] = df_bb1_gsp["GSP"].replace(gsp_mapping)
    extra_gsp_coordinates_df = (
        pd.DataFrame.from_dict(extra_gsp_coordinates, orient="index")
        .rename_axis(index="Name")
        .reset_index()
    )
    df_gsp_coordinates = (
        # There are cases of duplicate GSPs where the lat and lon information is the same but the GSP ID and GSP group are slightly different
        pd.concat([df_gsp_coordinates, extra_gsp_coordinates_df])
        .drop_duplicates(subset=["Name", "Latitude", "Longitude"])
        .dropna(subset=["Latitude", "Longitude"])
        .set_index("Name")
        .reset_index()
    )

    # Join GSP shape data with GSP coordinate data
    gsp_joined = df_gsp_shapes.set_index("GSPs").join(
        gdf_gsps.set_index("GSP ID"), lsuffix="shape", rsuffix="coord", how="outer"
    )

    # Calculate representative point and centroid for each GSP shape
    gsp_joined.loc[:, "representative_point"] = gsp_joined.geometryshape.to_crs(
        gdf_gsps.estimate_utm_crs()
    ).representative_point()
    gsp_joined.loc[:, "centroid"] = gsp_joined.geometryshape.to_crs(
        gdf_gsps.estimate_utm_crs()
    ).centroid

    fes_merged = pd.merge(
        df_bb1_gsp,
        gsp_joined.reset_index(),
        left_on="GSP",
        right_on="Name",
        how="outer",
    )

    # Drop unrequired columns
    fes_merged.drop(
        columns=["GSP Group", "Minor FLOP", "Latitude", "Longitude"], inplace=True
    )

    # Some GSPs have missing lat/lon information in the GSP coordinate file 
    unmatched = fes_merged.loc[
        (fes_merged.geometrycoord.isna()) & (fes_merged.GSPs.notna())
    ]
    for gsps in unmatched["GSPs"].tolist():
        all_occurences = fes_merged.loc[
            (fes_merged["GSPs"].str.contains(gsps)) & (fes_merged.GSPs.notna())
        ]
        if len(all_occurences) > 1:
            fes_merged = _merge_gsps(fes_merged, all_occurences, gsps)

    # Some busbars are split into multiple GSPs in the FES workbook but respresented as a single GSP in shape data
    for key in combine_busbars.keys():
        combine_busbars[key].append(key)
        gsps = "|".join(combine_busbars[key])
        all_occurences = fes_merged.loc[
            (fes_merged["GSPs"].str.contains(gsps)) & (fes_merged.GSPs.notna())
        ]
        if len(all_occurences) > 1:
            fes_merged = _merge_gsps(fes_merged, all_occurences, gsps)

    fes_merged = gpd.GeoDataFrame(fes_merged, geometry="geometryshape", crs="EPSG:4326")
    fes_merged[["representative_point", "centroid", "geometrycoord"]] = fes_merged[
        ["representative_point", "centroid", "geometrycoord"]
    ].to_wkt()
    return fes_merged


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the file paths
    bb1_path = snakemake.input.bb1_sheet
    bb2_path = snakemake.input.bb2_sheet
    gsp_coordinates_path = snakemake.input.gsp_coordinates
    gdf_regions = gpd.read_file(snakemake.input.regions)

    # Load all the params
    fes_scenario = get_scenario_name(snakemake)
    year_range = snakemake.params.year_range
    extra_gsp_coordinates = snakemake.params.fill_gsp_lat_lons
    manual_gsp_mapping = snakemake.params.manual_gsp_mapping
    combine_busbars = snakemake.params.combine_busbars

    df, df_bb1 = parse_inputs(
        bb1_path,
        bb2_path,
        gsp_coordinates_path,
        extra_gsp_coordinates,
        manual_gsp_mapping,
        fes_scenario,
        year_range,
    )

    shape = create_gsp_shapefile(
        gsp_coordinates_path,
        snakemake.input.gsp_shapes,
        df_bb1,
        gdf_regions,
        manual_gsp_mapping,
        extra_gsp_coordinates,
        combine_busbars,
    )

    region_data = map_points_to_regions(
        df,
        gdf_regions,
        "Latitude",
        "Longitude",
        "EPSG:4326",
        snakemake.params.target_crs,
    )[["name", "TO_region"]]
    df_with_regions = pd.concat(
        [df, region_data.rename(columns={"name": "bus"})], axis=1
    )
    for TO_region in gdf_regions["TO_region"].unique():
        df_with_regions.loc[
            df_with_regions.GSP == f"Direct({TO_region})", "TO_region"
        ] = TO_region
    if (null_bus := df_with_regions.bus.isnull()).any():
        warning_data = df_with_regions[null_bus][
            ["GSP", "Latitude", "Longitude", "TO_region"]
        ].drop_duplicates()
        logger.warning(
            f"There are GSPs with missing bus/region information after mapping lat/lon to regions:\n{warning_data}"
        )
    logger.info(f"Extracted the {fes_scenario} relevant data")

    df_with_regions.to_csv(snakemake.output.csv, index=False)

    shape.to_file(snakemake.output.shapefile, driver="GeoJSON")
