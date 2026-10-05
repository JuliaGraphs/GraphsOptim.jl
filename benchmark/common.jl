# Tools shared by all benchmark suites.

# All libraries use a single BLAS thread for a fair comparison.
# This must be set before Python and NumPy are loaded.
ENV["OPENBLAS_NUM_THREADS"] = "1"

using Graphs
using GraphsOptim
using HiGHS
using JuMP
using LinearAlgebra
using PythonCall
using Random
using SparseArrays
using Statistics

BLAS.set_num_threads(1)

const RESULTS_DIR = joinpath(@__DIR__, "results")
const FIGURES_DIR = joinpath(dirname(@__DIR__), "docs", "src", "assets", "benchmarks")

"Time limit (in seconds) given to the optimizers and to the exponential-time algorithms."
const TIME_LIMIT = 60.0

"HiGHS with the benchmark time limit."
const OPTIMIZER = optimizer_with_attributes(HiGHS.Optimizer, "time_limit" => TIME_LIMIT)

## Timing

"""
    measure(f; min_time=1.0, max_reps=7)

Run `f()` once (after compilation by a previous call) and repeat it until `min_time` seconds or `max_reps` runs are reached.
Return the median running time and the result of the last run, or `(NaN, nothing)` if `f` fails, e.g. because a time limit was hit.
"""
function measure(f; min_time=1.0, max_reps=7)
    times = Float64[]
    local result
    try
        while length(times) < max_reps && (isempty(times) || sum(times) < min_time)
            t = @elapsed result = f()
            push!(times, t)
        end
    catch e
        e isa InterruptException && rethrow()
        @warn "Run failed" exception = (e, catch_backtrace())
        return NaN, nothing
    end
    return median(times), result
end

pyexec(
    """
    import time

    def measure(f, args, min_time=1.0, max_reps=7, **kwargs):
        times = []
        while len(times) < max_reps and (not times or sum(times) < min_time):
            t = time.perf_counter()
            result = f(*args, **kwargs)
            times.append(time.perf_counter() - t)
        times.sort()
        n = len(times)
        median = times[n // 2] if n % 2 else (times[n // 2 - 1] + times[n // 2]) / 2
        return median, result
    """,
    Main,
)

"""
    pymeasure(f, args...; kwargs...)

Same as [`measure`](@ref) for the Python function `f`, timed on the Python side so that the conversion of the arguments is not counted.
"""
function pymeasure(f, args...; kwargs...)
    t, result = pyeval("measure", Main)(f, pytuple(args); kwargs...)
    return pyconvert(Float64, t), result
end

## Conversions

"Convert an undirected or directed `Graphs.jl` graph into a NetworkX graph with vertices `0:n-1`."
function to_networkx(g::AbstractGraph; vertex_weights=nothing)
    nx = pyimport("networkx")
    h = is_directed(g) ? nx.DiGraph() : nx.Graph()
    h.add_nodes_from(0:(nv(g) - 1))
    h.add_edges_from([(src(e) - 1, dst(e) - 1) for e in edges(g)])
    if !isnothing(vertex_weights)
        for v in vertices(g)
            h.nodes[v - 1]["weight"] = vertex_weights[v]
        end
    end
    return h
end

## Results

function write_csv(file, rows)
    mkpath(RESULTS_DIR)
    open(joinpath(RESULTS_DIR, file), "w") do io
        println(io, join(keys(first(rows)), ","))
        for row in rows
            println(io, join(values(row), ","))
        end
    end
    return nothing
end

function read_csv(file)
    lines = readlines(joinpath(RESULTS_DIR, file))
    header = Tuple(Symbol.(split(lines[1], ",")))
    parsevalue(x) = something(tryparse(Int, x), tryparse(Float64, x), String(x))
    return [
        NamedTuple{header}(Tuple(parsevalue.(split(line, ",")))) for line in lines[2:end]
    ]
end

function write_environment()
    commit = try
        readchomp(`git -C $(@__DIR__) describe --always --dirty`)
    catch
        "unknown"
    end
    version(mod) = pyconvert(String, pyimport(mod).__version__)
    open(joinpath(RESULTS_DIR, "environment.md"), "w") do io
        println(io, "- CPU: ", Sys.cpu_info()[1].model, " (single BLAS thread)")
        println(io, "- Julia ", VERSION, ", GraphsOptim.jl at commit `", commit, "`")
        println(
            io,
            "- Graphs.jl ",
            pkgversion(Graphs),
            ", JuMP.jl ",
            pkgversion(JuMP),
            ", HiGHS.jl ",
            pkgversion(HiGHS),
        )
        println(
            io,
            "- Python ",
            pyconvert(String, pyimport("platform").python_version()),
            ", graspologic ",
            version("graspologic"),
            ", NetworkX ",
            version("networkx"),
            ", SciPy ",
            version("scipy"),
        )
    end
    return nothing
end

"""
    too_slow(time)

Whether a method should be skipped on larger instances, because it failed or took more than a tenth of the time limit.
Exponential-time algorithms written in Python cannot be interrupted, so they are stopped early.
"""
too_slow(time) = isnan(time) || time > TIME_LIMIT / 10
