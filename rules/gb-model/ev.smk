# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Electric vehicle demand and V1G / V2G rules.
"""


rule process_ev_demand_shape:
    message:
        "Process EV demand profile shape into CSV format"
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        plug_in_offset=config["ev"]["plug_in_offset"],
        charging_duration=config["ev"]["charging_duration"],
    input:
        clustered_pop_layout=resources("pop_layout_base_s_clustered.csv"),
        traffic_data_KFZ="data/bundle/emobility/KFZ__count",
    output:
        demand_shape=resources("gb-model/ev_demand_shape.csv"),
    log:
        logs("process_ev_demand_shape.log"),
    script:
        "../../scripts/gb_model/ev/process_ev_demand_shape.py"


rule create_ev_peak_charging_table:
    message:
        "Process EV unmanaged charging demand from FES workbook into CSV format"
    params:
        scenario=config["fes"]["scenario"],
        year_range=config["fes"]["year_range_incl"],
    input:
        unmanaged_charging_sheet=resources("gb-model/fes/ED5.csv"),
    output:
        csv=resources("gb-model/ev_peak_demand.csv"),
    log:
        logs("create_ev_peak_charging_table.log"),
    script:
        "../../scripts/gb_model/ev/create_ev_peak_charging_table.py"


rule create_ev_v2g_storage_table:
    message:
        "Process EV unmanaged charging demand from FES workbook into CSV format"
    params:
        v2g_multiplier=config["fes"]["gb"]["flexibility"][
            "v2g_storage_to_capacity_ratio"
        ],
    input:
        v2g_cap=resources("gb-model/regional_ev_v2g.csv"),
    output:
        csv=resources("gb-model/regional_ev_v2g_storage.csv"),
    log:
        logs("create_ev_v2g_storage_table.log"),
    script:
        "../../scripts/gb_model/ev/create_ev_v2g_storage_table.py"


use rule scaled_demand_profile as scaled_ev_demand_profile with:
    input:
        gb_demand_annual=resources("gb-model/regional_{demand_type}_demand_annual.csv"),
        eur_demand_annual=resources("gb-model/eur_demand_annual.csv"),
        demand_shape=resources("gb-model/{demand_type}_demand_shape.csv"),
        gb_demand_peak=resources("gb-model/regional_{demand_type}_peak.csv"),
    params:
        scaling_params=config["ev"]["ev_demand_profile_transformation"],
    wildcard_constraints:
        demand_type="ev",
