# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT

import logging
import re

import numpy as np
import pandas as pd
import pypsa

logger = logging.getLogger(__name__)

FP = "../../data/fbmc/flow_based_constraints_v20260202.xlsx"

def load_fb_data(
    fp=FP # placeholder value
) -> pd.DataFrame:
    """
    Load PTDF/RAM matrix from Excel file.

    PTDF matrix contains the weights for each flow through each
    line/link, by each critical network element component (CNEC)
    known as a boundary in the GB PTDF.
    """

    # PTDF data is combined with RAM, etc.
    combined_data = pd.read_excel(FP)
    combined_data.columns = combined_data.columns.str.replace('^ptdf_', '', regex=True)
    combined_data = combined_data.melt(id_vars=['datetime', 'boundary name', 'direction', 'fref', 'f0', 'ram'], var_name='Link name', value_name='PTDF')

    # # unnecessary?
    # # Currently requires an existing reference network to map the buses using the NSIDE names. 
    # ref_n = pypsa.Network("~/Documents/oet/NGV-FBMC/modules/gb-dispatch-model/resources/GB/networks/unconstrained_clustered/2040.nc")
    # ref_n_links = ref_n.links[ref_n.links.carrier=='DC'][['bus0', 'bus1']]
    # combined_data = pd.merge(ref_n_links, combined_data, how='right', left_on='name', right_on='Link name') # do we need to keep the relation/ links?

    # # Flip the bus order so flow direction is positive
    # mask = ptdf_gb.direction == 'OPPOSITE'
    # combined_data.loc[mask, ['bus0', 'bus1']] = combined_data.loc[mask, ['bus1', 'bus0']].values

    # Rename/reindex the columns to match the existing FBMC 
    combined_data = combined_data.rename({'datetime':'snapshot', 'boundary name':'CNEC_ID','Link name':'name'}, axis=1)
    # johannes wants [OPP] and [DIR] in square brackets
    combined_data["direction"] = combined_data["direction"].apply(lambda x: f"[{x[:3]}]")
    combined_data["name"] = combined_data["name"] + " " + combined_data["direction"] # are opposing links handled in the og fbmc implementation?
    combined_data = combined_data.set_index(["CNEC_ID", "snapshot", "name"])

    return combined_data

def add_fbmc_constraints(n: pypsa.Network) -> None:
    """
    Add the FBMC constraints to the pypsa.Network model.

    Function is currently tailored towards the PTDF matrix and RAM values from ERAA2023,
    can be downloaded from https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/2023/FB-Domain-CORE_Merged.xlsx .

    Parameters
    ----------
    n : pypsa.Network
        The pypsa.Network object to which the FBMC constraints will be added.
    fp : str, optional
        File path to the Excel file containing the FBMC data.
        Needs to contain the PTDF matrix, RAM matrix, and weather assignments.
    config : dict
        Configuration used to modify the network for FBMC implementation.
    """

    fbmc_data = load_fb_data()

    # ----------------------------------
    # First part of the FBMC constraint:
    # Flows into and out of CORE bidding zones
    # ----------------------------------

    # get data for the individual study zone
    # fbmc_data[fbmc_data.index.get_level_values('name') == 'gb']

    fbmc_data_sz = fbmc_data[fbmc_data.index.get_level_values('name').isin(['gb [DIR]', 'gb [OPP]'])]

    # why isn't this working in the modify function?
    mask = n.components.links.static.bus1 == "EUR"
    n.components.links.static.loc[mask, "PTDF_type"] = "PTDF_SZ"

    # get links relevant for the intra-CCR FBMC constraint
    # this is GB # to GB CORE
    # links_idx = n.components.links.static.loc[
    #     n.components.links.static["PTDF_type"] == "PTDF_SZ"
    # ].index
    links_idx = n.components.links.static.loc[
        n.components.links.static["bus1"] == "EUR"
    ].index

    # # "study_zones" are named after the bidding zones, the flows in the network
    # # have a `-<CCR>` suffix, e.g. `-CORE` to indicate they are flows between the bidding zone
    # # and the CCR hub. Therefore, to align the PTDF data with the flows, we need to add the
    # # suffix to the study_zone names. We get the suffix from the link names and rename
    # # the study zones in the PTDF rather than the linopy model
    # rename_study_zones = {idx.split("-")[0]: idx for idx in links_idx}
    # if len(set(rename_study_zones.values())) != len(links_idx):
    #     raise ValueError(
    #         "Renaming of study zones to match link names resulted in a non 1:1 mapping."
    #     )
    # ptdf["study_zone"] = ptdf["study_zone"].replace(rename_study_zones)

    # TODO we do not include the offshore regions into the calculation
    # reason: Not part of the market, not interconnection to other countries,
    # consider with their transfer capacitiy rather than including them into the FBMC

    # get flow through links in CORE bidding zones

    # do the fancy multiplication
    ds = fbmc_data_sz["PTDF"].to_xarray()
    # Flows from GB # to virtual GB CORE (I think these flows need to be summed?)
    flows = n.model["Link-p"].sel(name=links_idx)
    mask = flows.coords["name"].str.endswith("ramp down")
    opp_flows = flows.sel(name=mask)    
    mask = flows.coords["name"].str.endswith("ramp up")
    dir_flows = flows.sel(name=mask)
    # third case where neither ramp down or ramp up is included?
    mask = ~flows.coords["name"].str.contains("ramp")
    net_flows = flows.sel(name=mask) + dir_flows

    # Calculate PTDF contribution and group by snapshot and CNEC_ID to sum up all contributions to each CNEC at each snapshot
    # need to think if this is the right orientation
    lhs_1_dir = ds.sel(name='gb [DIR]') * net_flows.sum(dim="name")
    lhs_1_opp = ds.sel(name='gb [OPP]') * opp_flows.sum(dim="name")
    lhs_1 = lhs_1_dir + lhs_1_opp
    # -----------------------------------
    # Second part of the FBMC constraint:
    # loading from HVDC lines between CORE and outside of CORE
    # -----------------------------------
    fbmc_data_ahc = fbmc_data.query('~name.str.contains("gb")')

    # Map pypsa.Network links that are related to DC and their names (index) to PTDF line names where bus0=from and bus1=to
    links = (
        n.components.links.static.query("`carrier`.str.startswith('DC')")[
            ["bus0", "bus1"]
        ]
        .reset_index()
        .rename(columns={"name": "link_name"})
    )

    ds = fbmc_data_ahc.to_xarray()
    mask = ds.coords["name"].str.endswith("[DIR]")
    net_ds = ds.sel(name=mask)
    mask = ds.coords["name"].str.endswith("[OPP]")
    opp_ds = ds.sel(name=mask)
    # Casting to xarray creates NaN values, need to fill those entries with 0
    # flows = n.model["Link-p"].sel(name=ds["name"]) # flows don't use the [OPP/DIR] syntax - they use ramp up/down
    # restructure flows
    breakpoint()
    net_ds = net_ds.reindex(name=net_flows.coords["name"])
    opp_ds = opp_ds.reindex(name=opp_flows.coords["name"])
    lhs_2_dir = net_ds * net_flows
    lhs_2_opp = opp_ds * opp_flows
    # Group by snapshot and CNEC_ID to sum up all contributions to each CNEC at each snapshot
    lhs_2 = (lhs_2_dir + lhs_2_opp).sum(dim="name")

    # # -----------------------------------
    # # Third part of the FBMC constraint:
    # # loading from HVDC lines within CORE region bidding zones
    # # -----------------------------------

    # # Load PTDF
    # ptdf = load_ptdf(fp=fp, sheet_name=ptdf_sheet_name, ptdf_type="PTDF_EvFB")

    # # Map PTDF to seasonal values for RAM
    # ptdf_snapshoted = wa.merge(
    #     ptdf,
    #     on="FB Domain",
    #     how="left",
    # )

    # flows = n.model["Link-p"].sel(name="EvFBA1-EvFBA2")

    # ds = ptdf_snapshoted.set_index(["CNEC_ID", "snapshot"])["PTDF"].to_xarray() * flows

    # ds = ds.reindex(
    #     snapshot=flows.coords["snapshot"],
    # )

    # # Casting to xarray creates NaN values, need to fill those entries with 0
    # lhs_3 = ds
    breakpoint()

    rhs = fbmc_data['ram'].droplevel("name").to_xarray()

    # Enable lhs_1 and lhs_3 when implemented
    n.model.add_constraints(
        lhs_1 + lhs_2 <= rhs,
        name="PTDF-RAM-constraints",
    )


def modify_network_for_fbmc(n: pypsa.Network) -> pypsa.Network:
    """
    Modify the pypsa.Network for the FBMC implementation.

    The methodology follows the description in ERAA2023.
    This function modified the network and adds additional components that are necessary for the
    evolved FBMC implementation.
    It also assigns some helpful, additional attributes to existing components like buses and links.

    Parameters
    ----------
    n : pypsa.Network
        The pypsa.Network object to be modified for FBMC implementation.
    config : dict
        Configuration used to modify the network for FBMC implementation.
        (TODO: Currently not used, but needed for proper ERAA2024 implementation)

    Returns
    -------
    pypsa.Network
        The modified pypsa.Network object with FBMC implementation.
    """

    # ---------------------------------------------------
    # Add the buses and links required for the FBMC
    # ---------------------------------------------------

    # 1. Virtual core

    # Any of these links is removed and replaced with an unlimited link
    # between the study zone and a virtual CORE hub
    n.add("Carrier", name="FBMC")
    n.add(
        "Bus",
        name="CORE GB",
        carrier="FBMC"
    )
    # links between the virtual GB specific core and the other GB nodes
    gb_buses = n.buses[n.buses.carrier == 'AC'].filter(like='GB', axis=0).index
    
    n.add(
        "Link",
        name=gb_buses,
        suffix='-GB CORE',
        bus0=gb_buses,
        bus1="CORE GB",
        p_nom=np.inf,
        efficiency=1.0,
        p_nom_extendable=False,
        p_min_pu=-1.0,
        p_max_pu=1.0,
        PTDF_type="PTDF_SZ", # not working
        FBMC_region="CORE",
    )

    mask = n.components.links.static.bus1 == "CORE GB"
    n.components.links.static.loc[mask, "PTDF_type"] = "PTDF_SZ"

    # ----------------------------------------------------
    # Add details on which FBMC region each bus belongs to
    # ----------------------------------------------------

    # Assign links an attribute to indicate which parts of the PTDF they are relevant for
    logger.info("Assigning PTDF types to network links for FBMC implementation.")
    n.components.buses.static["FBMC_region"] = ""
    n.components.links.static["PTDF_type"] = ""
    
    # 1. PTDF_SZ for intra-CORE flows
    core_buses = n.components.buses.static.query("FBMC_region == 'CORE'").index.tolist()
    idx = n.components.links.static[
        (n.components.links.static["bus0"].isin(core_buses))
        & (n.components.links.static["bus1"].isin(core_buses))
    ].index
    n.links.loc[idx, "PTDF_type"] = "PTDF_SZ"
    n.links.loc[idx, "FBMC_region"] = "CORE-CORE"

    # 2. PTDF*_AHC,SZ for flows between CORE and outside of CORE
    idx = n.components.links.static[
        (
            (n.components.links.static["bus0"].isin(core_buses))
            ^ (n.components.links.static["bus1"].isin(core_buses))
        )
        & (n.components.links.static["carrier"].isin(["DC", "DC_OH", "AC"]))
        & (n.components.links.static["PTDF_type"] != "PTDF_SZ")
    ].index
    n.links.loc[idx, "PTDF_type"] = "PTDF*_AHC,SZ"
    n.links.loc[idx, "FBMC_region"] = "CORE-Outside"

    return n