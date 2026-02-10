# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Building heat demand and demand-side response rules.
"""


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
        scripts("gb_model/heat/process_cop_profiles.py")


rule process_fes_heat_technologies:
    message:
        "Process the share of electrified heating technologies from FES workbook"
    params:
        year_range=config["redispatch"]["year_range_incl"],
        electrified_heating_technologies=config["fes"]["gb"]["demand"]["heat"][
            "electrified_heating_technologies"
        ],
        scenario=config["fes"]["scenario"],
    input:
        fes_heat_technology_data=resources(f"gb-model/fes/ED3.csv"),
    output:
        residential=resources("gb-model/residential_heat_techs_consumption.csv"),
        services=resources("gb-model/iandc_heat_techs_consumption.csv"),
    log:
        logs("process_fes_heat_technologies.log"),
    script:
        scripts("gb_model/heat/process_fes_heat_technologies.py")


rule resistive_heater_demand_profile:
    message:
        "Process resistive heat demand profile shape into CSV format for {wildcards.year}"
    input:
        energy_totals=resources("pop_weighted_energy_totals_s_clustered.csv"),
        heat_demand_shape=resources("hourly_heat_demand_total_base_s_clustered.nc"),
        residential_heat_techs_consumption=resources(
            "gb-model/regional_residential_heat_techs_consumption_inc_eur.csv"
        ),
        services_heat_techs_consumption=resources(
            "gb-model/regional_iandc_heat_techs_consumption_inc_eur.csv"
        ),
    output:
        csv=resources("gb-model/resistive_heater_demand/{year}.csv"),
    log:
        logs("resistive_heater_demand_profile_{year}.log"),
    script:
        scripts("gb_model/heat/resistive_heater_demand_profile.py")


rule process_heat_demand_shape:
    message:
        "Cluster default PyPSA-Eur {wildcards.sector} heat demand shape by bus for future year {wildcards.year}"
    input:
        demand=resources("hourly_heat_demand_total_base_s_clustered.nc"),
        cop_profile=resources("gb-model/cop/{year}.csv"),
        heating_mix=resources("gb-model/{sector}_heat_techs_consumption.csv"),
        energy_totals=resources("pop_weighted_energy_totals_s_clustered.csv"),
    output:
        csv=resources("gb-model/{sector}_heat_demand_shape/{year}.csv"),
    log:
        logs("heat_demand_s_clustered_{sector}_{year}.log"),
    script:
        scripts("gb_model/heat/process_heat_demand_shape.py")
