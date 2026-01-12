# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Hydrogen subsystem (demand, storage, electrolysis, fuel cells / turbine) rules.
"""


rule create_hydrogen_demand_table:
    message:
        "Process hydrogen demand data from FES workbook into CSV format"
    params:
        scenario=config["fes"]["gb"]["scenario"],
        year_range=config["fes"]["year_range_incl"],
        fes_demand_sheets=config["fes"]["hydrogen"]["demand"]["annual_demand_sheets"],
        other_sectors_list=config["fes"]["hydrogen"]["demand"]["other_sectors_list"],
    input:
        demand_sheets=lambda wildcards: [
            (
                resources(f"gb-model/fes/{config['fes-year']}/{sheet}.csv")
                if sheet != "other" or config["fes-year"] >= 2023
                else resources(f"gb-model/fes/2023/{sheet}.csv")
            )
            for sheet in config["fes"]["hydrogen"]["demand"]["annual_demand_sheets"][
                config["fes-year"]
            ]
        ],
    output:
        hydrogen_demand=resources("gb-model/fes_hydrogen_demand.csv"),
    log:
        logs("create_hydrogen_demand_table.log"),
    script:
        "../../scripts/gb_model/hydrogen/create_hydrogen_demand_table.py"


rule create_grid_electrolysis_table:
    message:
        "Process hydrogen electrolysis data from FES workbook into CSV format"
    input:
        regional_gb_data=resources("gb-model/regional_gb_data.csv"),
    output:
        csv=resources("gb-model/regional_grid_electrolysis_capacities.csv"),
    log:
        logs("create_grid_electrolysis_table.log"),
    script:
        "../../scripts/gb_model/hydrogen/create_grid_electrolysis_table.py"


rule create_hydrogen_supply_table:
    message:
        "Process hydrogen supply data from FES workbook into CSV format"
    params:
        scenario=config["fes"]["gb"]["scenario"],
        year_range=config["fes"]["year_range_incl"],
        fes_supply_sheets=config["fes"]["hydrogen"]["supply"]["supply_sheets"],
        exogeneous_supply_list=config["fes"]["hydrogen"]["supply"][
            "exogeneous_supply_list"
        ],
    input:
        supply_sheets=lambda wildcards: [
            resources(f"gb-model/fes/{year}/{sheet}.csv")
            for year, sheets in config["fes"]["hydrogen"]["supply"][
                "supply_sheets"
            ].items()
            for sheet in sheets.values()
        ],
    output:
        hydrogen_supply=resources("gb-model/fes_hydrogen_supply.csv"),
    log:
        logs("create_hydrogen_supply_table.log"),
    script:
        "../../scripts/gb_model/hydrogen/create_hydrogen_supply_table.py"


rule create_off_grid_electrolysis_demand:
    message:
        "Process electricity demand of off-grid electrolysis from FES workbook into CSV format"
    params:
        scenario=config["fes"]["gb"]["scenario"],
        year_range=config["fes"]["year_range_incl"],
        fes_supply_sheets=config["fes"]["hydrogen"]["supply"]["supply_sheets"],
    input:
        supply_sheets=lambda wildcards: [
            resources(f"gb-model/fes/{year}/{sheet}.csv")
            for year, sheets in config["fes"]["hydrogen"]["supply"][
                "supply_sheets"
            ].items()
            for sheet in sheets.values()
        ],
        grid_electrolysis_capacities=resources(
            "gb-model/regional_grid_electrolysis_capacities.csv"
        ),
    output:
        electricity_demand=resources(
            "gb-model/regional_off_grid_electrolysis_electricity_demand.csv"
        ),
    log:
        logs("create_off_grid_electrolysis_demand.log"),
    script:
        "../../scripts/gb_model/hydrogen/create_off_grid_electrolysis_demand.py"


rule create_hydrogen_storage_table:
    message:
        "Process hydrogen storage data from FES workbook into CSV format"
    params:
        scenario=config["fes"]["gb"]["scenario"],
        year_range=config["fes"]["year_range_incl"],
        fes_storage_sheets=config["fes"]["hydrogen"]["storage"]["storage_sheets"],
        interpolation_method=config["fes"]["hydrogen"]["storage"][
            "interpolation_method"
        ],
    input:
        storage_sheet=lambda wildcards: [
            resources(f"gb-model/fes/{year}/{sheet}.csv")
            for year, sheets in config["fes"]["hydrogen"]["storage"][
                "storage_sheets"
            ].items()
            for sheet in sheets.values()
        ],
    output:
        csv=resources("gb-model/H2_storage_capacity.csv"),
    log:
        logs("create_hydrogen_storage_table.log"),
    script:
        "../../scripts/gb_model/hydrogen/create_hydrogen_storage_table.py"


rule process_H2_demand:
    message:
        "Get net H2 demand"
    input:
        demand=resources("gb-model/fes_hydrogen_demand.csv"),
        fixed_supply=resources("gb-model/fes_hydrogen_supply.csv"),
    output:
        csv=resources("gb-model/H2_demand_annual.csv"),
    log:
        logs("process_H2_demand.log"),
    script:
        "../../scripts/gb_model/hydrogen/process_H2_demand.py"


rule process_regional_H2_data:
    message:
        "Process {wildcards.data_type} H2 data into GB regions format"
    input:
        to_regionalise=resources("gb-model/{data_type}.csv"),
        regional_distribution=resources(
            "gb-model/regional_grid_electrolysis_capacities.csv"
        ),
    params:
        param_name=lambda wildcards: {
            "H2_demand_annual": "p_set",
            "H2_storage_capacity": "e_nom",
        }[wildcards.data_type],
    output:
        csv=resources("gb-model/regional_{data_type}.csv"),
    wildcard_constraints:
        data_type="H2_demand_annual|H2_storage_capacity",
    log:
        logs("process_regional_H2_data_{data_type}.log"),
    script:
        "../../scripts/gb_model/hydrogen/process_regional_H2_data.py"


rule add_eur_H2_demand:
    message:
        "Add European H2 demand based on historical data combined with TYNDP future scenario demands"
    input:
        gb_demand=resources("gb-model/regional_H2_demand_annual.csv"),
        eur_demand_tyndp="data/gb-model/tyndp_h2_demand.csv",
        eur_demand_today="data/gb-model/downloaded/eur_H2_demand_today.xlsx",
    params:
        countries=config["countries"],
    output:
        csv=resources("gb-model/regional_H2_demand_annual_inc_eur.csv"),
    log:
        logs("add_eur_H2_demand.log"),
    script:
        "../../scripts/gb_model/hydrogen/add_eur_H2_demand.py"


rule synthesise_eur_H2_data:
    message:
        "Synthesise European H2 {wildcards.dataset} data using GB data"
    input:
        h2_demand=resources("gb-model/regional_H2_demand_annual_inc_eur.csv"),
        gb_only_dataset=resources("gb-model/regional_{dataset}.csv"),
    output:
        csv=resources("gb-model/regional_{dataset}_inc_eur.csv"),
    log:
        logs("synthesise_eur_H2_data_{dataset}.log"),
    wildcard_constraints:
        dataset="off_grid_electrolysis_electricity_demand|H2_storage_capacity|grid_electrolysis_capacities",
    script:
        "../../scripts/gb_model/hydrogen/synthesise_eur_H2_data.py"
