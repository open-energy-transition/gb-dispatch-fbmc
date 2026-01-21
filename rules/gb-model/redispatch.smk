# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Boundary constrained redispatch run rules.
"""


rule process_CfD_strike_prices:
    message:
        "get strike price for low carbon contracts"
    params:
        carrier_mapping=config["low_carbon_register"]["carrier_mapping"],
        end_year=config["redispatch"]["year_range_incl"][0],
    input:
        register="data/gb-model/downloaded/low-carbon-contracts.csv",
    output:
        csv=resources("gb-model/CfD_strike_prices.csv"),
    log:
        logs("process_CfD_strike_prices.log"),
    script:
        "../../scripts/gb_model/redispatch/process_CfD_strike_prices.py"


rule calc_interconnector_bid_offer_profile:
    message:
        "Calculate interconnector bid/offer profiles"
    params:
        bids_and_offers=config_provider("redispatch"),
    input:
        unconstrained_result=RESULTS + "networks/unconstrained_clustered/{year}.nc",
    output:
        bid_offer_profile=resources(
            "gb-model/bids_and_offers/{year}/interconnector_bid_offer_profile.csv"
        ),
    log:
        logs("calc_interconnector_bid_offer_profile/{year}.log"),
    script:
        "../../scripts/gb_model/redispatch/calc_interconnector_bid_offer_profile.py"


rule prepare_constrained_network:
    message:
        "Prepare network for constrained optimization"
    params:
        bids_and_offers=config_provider("redispatch"),
    input:
        network=resources("networks/composed_clustered/{year}.nc"),
        unconstrained_result=RESULTS + "networks/unconstrained_clustered/{year}.nc",
        renewable_strike_prices=resources("gb-model/CfD_strike_prices.csv"),
        interconnector_bid_offer=resources(
            "gb-model/bids_and_offers/{year}/interconnector_bid_offer_profile.csv"
        ),
    output:
        network=resources("networks/constrained_clustered/{year}.nc"),
    log:
        logs("prepare_constrained_network/{year}.log"),
    script:
        "../../scripts/gb_model/redispatch/prepare_constrained_network.py"


rule solve_constrained:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=Path(workflow.snakefile).parent
        / "../../scripts/gb_model/redispatch/custom_constraints.py",
        etys_boundaries_to_lines=config["region_operations"]["etys_boundaries_lines"],
        etys_boundaries_to_links=config["region_operations"]["etys_boundaries_links"],
    input:
        network=resources("networks/constrained_clustered/{year}.nc"),
        etys_caps=resources("gb-model/etys_boundary_capabilities.csv"),
    output:
        network=RESULTS + "networks/constrained_clustered/{year}.nc",
        config=RESULTS + "configs/config.constrained_clustered/{year}.yaml",
    log:
        solver=normpath(
            RESULTS + "logs/solve_network/constrained_clustered/{year}_solver.log"
        ),
        memory=RESULTS + "logs/solve_network/constrained_clustered/{year}_memory.log",
        python=RESULTS + "logs/solve_network/constrained_clustered/{year}_python.log",
    benchmark:
        (RESULTS + "benchmarks/solve_network/constrained_clustered/{year}")
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="6h"),
    shadow:
        shadow_config
    script:
        "../../scripts/solve_network.py"


rule calculate_constraint_costs:
    message:
        "Calculate total constraint cost across all years"
    params:
        constraint_cost_extra_years=config["redispatch"]["constraint_cost_extra_years"],
    input:
        networks=expand(
            RESULTS + "networks/constrained_clustered/{year}.nc",
            year=range(
                config["redispatch"]["year_range_incl"][0],
                config["redispatch"]["year_range_incl"][1] + 1,
            ),
        ),
    output:
        csv=RESULTS + "constraint_cost.csv",
    script:
        "../../scripts/gb_model/redispatch/calculate_constraint_costs.py"
