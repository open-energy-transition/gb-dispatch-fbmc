# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Storage unit component rules.
"""


rule process_regional_battery_storage_capacity:
    message:
        "Process national storage data from FES workbook into CSV format"
    params:
        year_range=config["redispatch"]["year_range_incl"],
        carrier_mapping=config["fes"]["gb"]["flexibility"]["carrier_mapping"]["battery"],
    input:
        flexibility_sheet=resources(f"gb-model/fes/FLX1.csv"),
        regional_data=resources("gb-model/{fes_scenario}/fes_powerplants.csv"),
    output:
        csv=resources(
            "gb-model/{fes_scenario}/regional_battery_storage_capacity_inc_eur.csv"
        ),
    log:
        logs("process_regional_battery_storage_capacity_{fes_scenario}.log"),
    script:
        scripts("gb_model/storage/process_regional_battery_storage_capacity.py")
