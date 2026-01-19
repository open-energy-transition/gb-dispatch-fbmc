# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Unconstrained dispatch run rules.
"""



rule prepare_unconstrained:
    message:
        "Prepare network for unconstrained optimization"
    input:
        network=resources("networks/composed_clustered/{year}.nc"),
    output:
        network=resources("networks/unconstrained_clustered/{year}.nc"),
    params:
        load_shedding_cost_above_marginal=config["fes"]["eur"][
            "load_shedding_cost_above_marginal"
        ],
    log:
        logs("prepare_unconstrained_network/{year}.log"),
    script:
        "../../scripts/gb_model/dispatch/prepare_unconstrained_network.py"


rule solve_unconstrained:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=[],
    input:
        network=resources("networks/unconstrained_clustered/{year}.nc"),
    output:
        network=RESULTS + "networks/unconstrained_clustered/{year}.nc",
        config=RESULTS + "configs/config.unconstrained_clustered/{year}.yaml",
    log:
        solver=normpath(
            RESULTS + "logs/solve_network/unconstrained_clustered/{year}_solver.log"
        ),
        memory=RESULTS + "logs/solve_network/unconstrained_clustered/{year}_memory.log",
        python=RESULTS + "logs/solve_network/unconstrained_clustered/{year}_python.log",
    benchmark:
        (RESULTS + "benchmarks/solve_network/unconstrained_clustered/{year}")
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="6h"),
    shadow:
        shadow_config
    script:
        "../../scripts/solve_network.py"
