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
        scenario=config["fes"]["gb"]["scenario"],
        year_range=config["fes"]["year_range_incl"],
        carrier_mapping=config["fes"]["gb"]["flexibility"]["carrier_mapping"]["battery"],
    input:
        flexibility_sheet=resources("gb-model/fes/2021/FLX1.csv"),
        regional_data=resources("gb-model/fes_powerplants.csv"),
    output:
        csv=resources("gb-model/regional_battery_storage_capacity_inc_eur.csv"),
    log:
        logs("process_regional_battery_storage_capacity.log"),
    script:
        "../../scripts/gb_model/storage/process_regional_battery_storage_capacity.py"
