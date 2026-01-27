# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Baseline and generalised demand and demand-side response (flexibility) rules.
"""


rule create_demand_table:
    message:
        "Process {wildcards.demand_type} demand from FES workbook into CSV format"
    params:
        technology_detail=config["fes"]["gb"]["demand"]["Technology Detail"],
    input:
        regional_gb_data=resources("gb-model/regional_gb_data.csv"),
    output:
        demand=resources("gb-model/regional_{demand_type}_demand_annual.csv"),
    wildcard_constraints:
        demand_type="|".join(config["fes"]["gb"]["demand"]["Technology Detail"]),
    log:
        logs("create_{demand_type}_demand_table.log"),
    script:
        "../../scripts/gb_model/demand_and_dsr/create_demand_table.py"


rule cluster_baseline_electricity_demand_timeseries:
    message:
        "Cluster default PyPSA-Eur baseline electricity demand timeseries by bus"
    params:
        scaling_factor=config_provider("load", "scaling_factor"),
    input:
        load=resources("electricity_demand_base_s.nc"),
        busmap=resources("busmap_base_s_clustered.csv"),
    output:
        csv_file=resources(
            "gb-model/historical_baseline_electricity_demand_profile.csv"
        ),
    log:
        logs("cluster_baseline_electricity_demand_timeseries.log"),
    script:
        "../../scripts/gb_model/demand_and_dsr/cluster_baseline_electricity_demand_timeseries.py"


rule process_baseline_demand_shape:
    message:
        "Process baseline electricity demand profile shape into CSV format for {wildcards.year}"
    input:
        historical_profile=resources(
            "gb-model/historical_baseline_electricity_demand_profile.csv"
        ),
        energy_totals=resources("pop_weighted_energy_totals_s_clustered.csv"),
        heat_demand_shape=resources("hourly_heat_demand_total_base_s_clustered.nc"),
        resistive_heater_demand=resources("gb-model/resistive_heater_demand/{year}.csv"),
    output:
        demand_shape=resources("gb-model/baseline_electricity_demand_shape/{year}.csv"),
    log:
        logs("process_baseline_demand_shape_{year}.log"),
    script:
        "../../scripts/gb_model/demand_and_dsr/process_baseline_demand_shape.py"


rule create_flexibility_table:
    message:
        "Process {wildcards.flexibility_type} flexibility from FES workbook into CSV format"
    params:
        scenario=config["fes"]["scenario"],
        year_range=config["redispatch"]["year_range_incl"],
        carrier_mapping=lambda wildcards: config["fes"]["gb"]["flexibility"][
            "carrier_mapping"
        ][wildcards.flexibility_type],
    input:
        flexibility_sheet=resources(f"gb-model/fes/FLX1.csv"),
    output:
        flexibility=resources("gb-model/{flexibility_type}.csv"),
    log:
        logs("create_flexibility_table_{flexibility_type}.log"),
    wildcard_constraints:
        flexibility_type="|".join(config["fes"]["gb"]["flexibility"]["carrier_mapping"]),
    script:
        "../../scripts/gb_model/demand_and_dsr/create_flexibility_table.py"


rule synthesise_gb_regional_data:
    message:
        "Process regional {wildcards.data} data from FES workbook into CSV format"
    params:
        regional_distribution_reference=lambda wildcards: config["fes"]["gb"][
            "regional_distribution_reference"
        ][wildcards.data],
    input:
        national_gb_data=resources("gb-model/{data}.csv"),
        regional_gb_data=resources("gb-model/regional_gb_data.csv"),
    output:
        csv=resources("gb-model/regional_{data}.csv"),
    wildcard_constraints:
        data="|".join(config["fes"]["gb"]["regional_distribution_reference"].keys()),
    log:
        logs("synthesise_gb_regional_data_{data}.log"),
    script:
        "../../scripts/gb_model/demand_and_dsr/synthesise_gb_regional_data.py"


rule distribute_eur_demands:
    message:
        "Distribute total European neighbour annual demands into base electricity, heating, and transport"
    input:
        eur_data=resources("gb-model/national_eur_data.csv"),
        energy_totals=resources("energy_totals.csv"),
        electricity_demands=expand(
            resources("gb-model/regional_{demand_type}_demand_annual.csv"),
            demand_type=config["fes"]["gb"]["demand"]["Technology Detail"].keys(),
        ),
        extra_demands=[],
    params:
        totals_to_demands=config["fes"]["eur"]["totals_to_demand_groups"],
        base_year=config["energy"]["energy_totals_year"],
    output:
        csv=resources("gb-model/eur_demand_annual.csv"),
    log:
        logs("distribute_eur_demands.log"),
    script:
        "../../scripts/gb_model/demand_and_dsr/distribute_eur_demands.py"


def _ref_demand_type(w):
    return config["fes"]["eur"]["add_data_reference"][w.dataset]


rule synthesise_eur_data:
    message:
        "Create a regional {wildcards.dataset} dataset including European neighbours based on GB data and relative annual demand"
    input:
        gb_demand_annual=lambda wildcards: resources(
            f"gb-model/regional_{_ref_demand_type(wildcards)}_demand_annual.csv"
        ),
        eur_demand_annual=resources("gb-model/eur_demand_annual.csv"),
        gb_only_dataset=resources("gb-model/regional_{dataset}.csv"),
    params:
        demand_type=_ref_demand_type,
    output:
        csv=resources("gb-model/regional_{dataset}_inc_eur.csv"),
    log:
        logs("synthesise_eur_data_{dataset}.log"),
    wildcard_constraints:
        dataset="|".join(config["fes"]["eur"]["add_data_reference"].keys()),
    script:
        "../../scripts/gb_model/demand_and_dsr/synthesise_eur_data.py"


rule scaled_demand_profile:
    message:
        "Generate {wildcards.demand_type} demand profile for model year {wildcards.year}"
    input:
        gb_demand_annual=resources("gb-model/regional_{demand_type}_demand_annual.csv"),
        eur_demand_annual=resources("gb-model/eur_demand_annual.csv"),
        demand_shape=resources("gb-model/{demand_type}_demand_shape/{year}.csv"),
    output:
        csv=resources("gb-model/{demand_type}_demand/{year}.csv"),
    log:
        logs("scaled_demand_profile_{demand_type}_{year}.log"),
    wildcard_constraints:
        demand_type="baseline_electricity|residential_heat|iandc_heat",
    script:
        "../../scripts/gb_model/demand_and_dsr/scaled_demand_profile.py"
