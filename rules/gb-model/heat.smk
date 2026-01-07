# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Building heat demand and demand-side response rules.
"""


rule create_heat_flexibility_table:
    message:
        "Process residential heat demand flexibility from FES workbook"
    params:
        scenario=config["fes"]["gb"]["scenario"],
        year_range=config["fes"]["year_range_incl"],
        flex_level=config["fes"]["gb"]["flexibility"]["residential_heat_flex_level"],
    input:
        flexibility_sheet=resources("gb-model/fes/2021/FL.10.csv"),
    output:
        csv=resources("gb-model/residential_heat_dsr_flexibility.csv"),
    log:
        logs("create_heat_flexibility_table.log"),
    script:
        "../../scripts/gb_model/heat/create_heat_flexibility_table.py"


rule process_cop_profiles:
    message:
        "Process COP profile for {wildcards.year} obtained from existing PyPSA-Eur rules"
    params:
        year=lambda wildcards: wildcards.year,
        heat_pump_sources=config["sector"]["heat_pump_sources"],
    input:
        cop_profile=resources("cop_profiles_base_s_clustered_{year}.nc"),
        clustered_pop_layout=resources("pop_layout_base_s_clustered.csv"),
        district_heat_share=resources("district_heat_share.csv"),
    output:
        csv=resources("gb-model/cop/{year}.csv"),
    log:
        logs("process_cop_profiles_clustered_{year}.log"),
    script:
        "../../scripts/gb_model/heat/process_cop_profiles.py"


rule process_fes_heating_mix:
    message:
        "Process the share of electrified heating technologies from FES workbook"
    params:
        year_range=config["fes"]["year_range_incl"],
        electrified_heating_technologies=config["fes"]["gb"]["demand"]["heat"][
            "electrified_heating_technologies"
        ],
        scenario=config["fes"]["gb"]["scenario"],
    input:
        fes_residential_heatmix=resources("gb-model/fes/2021/CV.16.csv"),
        fes_services_heatmix=resources("gb-model/fes/2021/CV.55.csv"),
        fes_hp_uptake_trend=resources("gb-model/fes/2021/CV.14.csv"),
    output:
        csv=resources("gb-model/fes_heating_mix.csv"),
    log:
        logs("process_fes_heating_mix.log"),
    script:
        "../../scripts/gb_model/heat/process_fes_heating_mix.py"


rule process_heat_demand_shape:
    message:
        "Cluster default PyPSA-Eur heat demand shape by bus"
    input:
        demand=resources("hourly_heat_demand_total_base_s_clustered.nc"),
        cop_profile=resources("gb-model/cop/{year}.csv"),
        heating_mix=resources("gb-model/fes_heating_mix.csv"),
        energy_totals=resources("pop_weighted_energy_totals_s_clustered.csv"),
    output:
        residential_csv=resources("gb-model/residential_heat_demand_shape/{year}.csv"),
        #Industry load is not generated in PyPSA-Eur, hence the same profile as services is considered to be applicable for c&i
        services_csv=resources("gb-model/iandc_heat_demand_shape/{year}.csv"),
    log:
        logs("heat_demand_s_clustered_{year}.log"),
    script:
        "../../scripts/gb_model/heat/process_heat_demand_shape.py"


rule resistive_heat_demand:
    message:
        "Get resistive electrical annual heat supply."
    input:
        gb_residential_heat_annual=resources(
            "gb-model/regional_residential_heat_demand_annual.csv"
        ),
        gb_iandc_heat_annual=resources("gb-model/regional_iandc_heat_demand_annual.csv"),
        eur_demand_annual=resources("gb-model/eur_demand_annual.csv"),
        heating_mix=resources("gb-model/fes_heating_mix.csv"),
    output:
        csv=resources("gb-model/resistive_heat_demand.csv"),
    log:
        logs("resistive_heat_demand.log"),
    script:
        "../../scripts/gb_model/heat/resistive_heat_demand.py"


use rule scaled_demand_profile as scaled_heat_demand_profile with:
    input:
        gb_demand_annual=resources("gb-model/{demand_type}_demand_annual.csv"),
        eur_demand_annual=resources("gb-model/eur_demand_annual.csv"),
        demand_shape=resources("gb-model/{demand_type}_demand_shape/{year}.csv"),
    wildcard_constraints:
        demand_type="residential_heat|iandc_heat",
