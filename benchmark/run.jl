# Benchmarks of GraphsOptim.jl against other Julia and Python libraries.
#
# Usage, from the root of the repository:
#
#     julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
#     julia --project=benchmark benchmark/run.jl [suite...] [--plots-only]
#
# Without arguments, all suites are run. The results are written to `benchmark/results/`
# and the figures to `docs/src/assets/benchmarks/`, both shown in the documentation.

include("common.jl")
include("plots.jl")
include("suites/graph_matching.jl")
include("suites/edit_distance.jl")
include("suites/paths_flows.jl")
include("suites/cliques.jl")

const SUITES = Dict(
    "graph_matching" => (run_graph_matching, plot_graph_matching),
    "edit_distance" => (run_edit_distance, plot_edit_distance),
    "paths_flows" => (run_paths_flows, plot_paths_flows),
    "cliques" => (run_cliques, plot_cliques),
)

function main(args)
    plots_only = "--plots-only" in args
    names = filter(!=("--plots-only"), args)
    isempty(names) && (names = sort(collect(keys(SUITES))))
    for name in names
        haskey(SUITES, name) || error(
            "Unknown suite $name, available suites: $(join(sort(collect(keys(SUITES))), ", "))",
        )
    end
    mkpath(FIGURES_DIR)
    for name in names
        run_suite, plot_suite = SUITES[name]
        if !plots_only
            @info "Running suite $name"
            run_suite()
        end
        plot_suite()
    end
    plots_only || write_environment()
    return nothing
end

main(ARGS)
