# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Generator component rules.
"""


rule retrieve_entsoe_unavailability_data:
    message:
        "Retrieve data from ENTSOE API for generator {wildcards.business_type} unavailability in {wildcards.zone} bidding zone"
    output:
        xml_base_dir=directory("data/gb-model/entsoe_api/{zone}/{business_type}"),
    params:
        start_date=config["entsoe_unavailability"]["start_date"],
        end_date=config["entsoe_unavailability"]["end_date"],
        bidding_zones=config["entsoe_unavailability"]["bidding_zones"],
        business_types=config["entsoe_unavailability"]["business_types"],
        max_request_days=config["entsoe_unavailability"]["max_request_days"],
        api_params=config["entsoe_unavailability"]["api_params"],
    log:
        logs("retrieve_entsoe_unavailability_data_{zone}_{business_type}.log"),
    resources:
        mem_mb=1000,
    script:
        "../../scripts/gb_model/generators/retrieve_entsoe_unavailability_data.py"


rule process_entsoe_unavailability_data:
    input:
        xml_base_dir="data/gb-model/entsoe_api/{zone}/{business_type}",
    output:
        unavailability=resources(
            "gb-model/{zone}_{business_type}_generator_unavailability.csv"
        ),
    log:
        logs("process_entsoe_unavailability_data_{zone}_{business_type}.log"),
    params:
        business_type_codes=config["entsoe_unavailability"]["api_params"][
            "business_types"
        ],
    resources:
        mem_mb=1000,
    script:
        "../../scripts/gb_model/generators/process_entsoe_unavailability_data.py"


rule generator_monthly_availability_fraction:
    input:
        planned=resources("gb-model/{zone}_planned_generator_unavailability.csv"),
        forced=resources("gb-model/{zone}_forced_generator_unavailability.csv"),
        powerplants=resources("powerplants_s_all.csv"),
    params:
        carrier_mapping=config["entsoe_unavailability"]["carrier_mapping"],
        resource_type_mapping=config["entsoe_unavailability"]["resource_type_mapping"],
        start_date=config["entsoe_unavailability"]["start_date"],
        end_date=config["entsoe_unavailability"]["end_date"],
        max_unavailable_days=config["entsoe_unavailability"]["max_unavailable_days"],
    output:
        csv=resources("gb-model/{zone}_generator_monthly_availability_fraction.csv"),
    log:
        logs("{zone}_generator_monthly_availability_fraction.log"),
    script:
        "../../scripts/gb_model/generators/generator_monthly_availability_fraction.py"


rule create_powerplants_table:
    message:
        "Tabulate powerplant data GSP-wise from FES workbook sheet BB1 and EU supply data"
    params:
        gb_config=config["fes"]["gb"],
        eur_config=config["fes"]["eur"],
        dukes_config=config["dukes-5.11"],
        default_set=config["fes"]["default_set"],
    input:
        gsp_data=resources("gb-model/regional_gb_data.csv"),
        eur_data=resources("gb-model/national_eur_data.csv"),
        dukes_data=resources("gb-model/dukes-current-capacity.csv"),
    output:
        csv=resources("gb-model/fes_powerplants.csv"),
    log:
        logs("create_powerplants_table.log"),
    script:
        "../../scripts/gb_model/generators/create_powerplants_table.py"


rule assign_costs:
    message:
        "Prepares costs file from technology-data of PyPSA-Eur and FES and assigns to {wildcards.data_file}"
    params:
        costs_config=config["costs"],
        fes_scenario=config["fes"]["scenario"],
    input:
        tech_costs=lambda w: resources(
            f"costs_{config_provider('costs', 'year')(w)}.csv"
        ),
        fes_power_costs=resources("gb-model/fes-costing/AS.1 (Power Gen).csv"),
        fes_carbon_costs=resources("gb-model/fes-costing/AS.7 (Carbon Cost).csv"),
        fes_powerplants=resources("gb-model/{data_file}.csv"),
    output:
        enriched_powerplants=resources("gb-model/{data_file}_inc_tech_data.csv"),
    log:
        logs("assign_costs_{data_file}.log"),
    wildcard_constraints:
        data_file="fes_powerplants|regional_H2_storage_capacity_inc_eur|regional_grid_electrolysis_capacities_inc_eur",
    script:
        "../../scripts/gb_model/generators/assign_costs.py"


rule create_chp_p_min_pu_profile:
    message:
        "Create CHP minimum operation profiles linked to heat demand"
    params:
        heat_to_power_ratio=config["chp"]["heat_to_power_ratio"],
        min_operation_level=config["chp"]["min_operation_level"],
        shutdown_threshold=config["chp"]["shutdown_threshold"],
    input:
        regions=resources("gb-model/merged_shapes.geojson"),
        heat_demand=resources("hourly_heat_demand_total_base_s_clustered.nc"),
    output:
        chp_p_min_pu=resources("gb-model/chp_p_min_pu.csv"),
    log:
        logs("create_chp_p_min_pu_profile.log"),
    script:
        "../../scripts/gb_model/generators/create_chp_p_min_pu_profile.py"
