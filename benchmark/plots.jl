# Figures of the benchmark results, saved as SVG files for the documentation.

using CairoMakie

const SURFACE = colorant"#fcfcfb"
const TEXT_PRIMARY = colorant"#0b0b0b"
const TEXT_SECONDARY = colorant"#52514e"
const GRID = colorant"#e4e3df"

# Categorical palette, in fixed order. The color of a series follows its library.
const BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
const MAGENTA, GREEN, VIOLET, RED = "#e87ba4", "#008300", "#4a3aa7", "#e34948"

"Color and line style of each series: heuristics share the color of their library, dashed."
const SERIES_STYLE = Dict(
    "GraphsOptim.jl" => (BLUE, :solid),
    "GraphsOptim.jl (native)" => (BLUE, :solid),
    "GraphsOptim.jl (JuMP)" => (VIOLET, :solid),
    "graspologic" => (ORANGE, :solid),
    "NetworkX" => (ORANGE, :solid),
    "NetworkX (approximation)" => (ORANGE, :dash),
    "NetworkX (2-approximation)" => (ORANGE, :dash),
    "Graphs.jl" => (AQUA, :solid),
    "Graphs.jl (all maximal cliques)" => (AQUA, :solid),
    "Graphs.jl (greedy)" => (AQUA, :dash),
    "Graphs.jl (2-approximation)" => (AQUA, :dot),
    "SciPy" => (YELLOW, :solid),
    # graph edit distance formulations
    "F1" => (BLUE, :solid),
    "F1prime" => (AQUA, :solid),
    "F1plus" => (YELLOW, :solid),
    "F2" => (MAGENTA, :solid),
    "F2minus" => (GREEN, :solid),
    "F2plus" => (VIOLET, :solid),
    "FORI" => (RED, :solid),
)

const THEME = Theme(;
    backgroundcolor=SURFACE,
    textcolor=TEXT_PRIMARY,
    fontsize=13,
    Axis=(;
        width=280,
        height=220,
        backgroundcolor=SURFACE,
        xgridcolor=GRID,
        ygridcolor=GRID,
        xgridwidth=1,
        ygridwidth=1,
        xminorgridvisible=false,
        yminorgridvisible=false,
        topspinevisible=false,
        rightspinevisible=false,
        leftspinecolor=TEXT_SECONDARY,
        bottomspinecolor=TEXT_SECONDARY,
        xtickcolor=TEXT_SECONDARY,
        ytickcolor=TEXT_SECONDARY,
        xticklabelcolor=TEXT_SECONDARY,
        yticklabelcolor=TEXT_SECONDARY,
        titlesize=14,
        titlefont=:bold,
    ),
    Legend=(; framevisible=false, labelcolor=TEXT_PRIMARY, backgroundcolor=SURFACE),
)

"""
    series!(ax, rows; x, y, group, order, halo)

Plot one line per series of `rows` (grouped by `group`), with `x` and `y` given as functions of a row.
The series listed in `halo` are drawn first as wide translucent bands, so that the other series stay visible when they coincide.
"""
function series!(ax, rows; x, y, group, order=unique(group.(rows)), halo=())
    for name in sort(order; by=name -> !(name in halo))
        rs = sort(filter(r -> group(r) == name, rows); by=x)
        xs, ys = Float64.(x.(rs)), Float64.(y.(rs))
        keep = .!isnan.(ys)
        any(keep) || continue
        color, linestyle = SERIES_STYLE[name]
        if name in halo
            lines!(ax, xs[keep], ys[keep]; label=name, color=(color, 0.35), linewidth=8)
        else
            scatterlines!(
                ax,
                xs[keep],
                ys[keep];
                label=name,
                color,
                linestyle,
                linewidth=2,
                markersize=8,
                strokecolor=SURFACE,
                strokewidth=2,
            )
        end
    end
    readable_ticks!(ax, rows; x, y)
    return nothing
end

function legend_below!(fig, ax; nbanks=1)
    return Legend(
        fig[end + 1, :],
        ax;
        orientation=:horizontal,
        nbanks,
        tellheight=true,
        tellwidth=false,
    )
end

function format_time(t)
    t >= 1 && return "$(round(Int, t)) s"
    t >= 1e-3 && return "$(round(Int, 1e3t)) ms"
    return "$(round(Int, 1e6t)) µs"
end

"Readable ticks on logarithmic axes: the actual sizes horizontally, decades of seconds vertically."
function readable_ticks!(ax, rows; x, y)
    if ax.xscale[] === log10
        xs = sort(unique(x.(rows)))
        ax.xticks = (Float64.(xs), string.(xs))
    end
    if ax.yscale[] === log10
        ys = filter(!isnan, Float64.(y.(rows)))
        decades = 10.0 .^ (floor(Int, log10(minimum(ys))):ceil(Int, log10(maximum(ys))))
        # with few decades, add intermediate ticks
        ticks = if maximum(ys) / minimum(ys) < 300
            [d * m for d in decades for m in (1, 2, 5)]
        else
            decades
        end
        ax.yticks = (ticks, format_time.(ticks))
    end
    return nothing
end

seconds_axis(pos; kwargs...) = Axis(pos; yscale=log10, ylabel="time", kwargs...)

function save_figure(name, fig)
    resize_to_layout!(fig)
    return save(joinpath(FIGURES_DIR, name), fig)
end

function plot_graph_matching()
    scaling = read_csv("graph_matching_scaling.csv")
    quality = read_csv("graph_matching_quality.csv")
    with_theme(THEME) do
        fig = Figure()
        for (k, algorithm) in enumerate(("FAQ", "GOAT"))
            ax = seconds_axis(
                fig[1, k]; title=algorithm, xscale=log10, xlabel="number of vertices n"
            )
            series!(
                ax,
                filter(r -> r.algorithm == algorithm, scaling);
                x=r -> r.n,
                y=r -> r.time,
                group=r -> r.library,
            )
            k == 1 && legend_below!(fig, ax)
        end
        save_figure("graph_matching_time.svg", fig)

        fig = Figure()
        for (k, algorithm) in enumerate(("FAQ", "GOAT"))
            ax = Axis(
                fig[1, k];
                title=algorithm,
                xlabel="correlation ρ",
                ylabel="accuracy",
                limits=(nothing, (-0.05, 1.05)),
            )
            rows = [
                (;
                    ρ,
                    library,
                    accuracy=mean(
                        r.accuracy for r in quality if
                        r.ρ == ρ && r.library == library && r.algorithm == algorithm
                    ),
                ) for ρ in unique(r.ρ for r in quality) for
                library in ("GraphsOptim.jl", "graspologic")
            ]
            series!(
                ax,
                rows;
                x=r -> r.ρ,
                y=r -> r.accuracy,
                group=r -> r.library,
                halo=("graspologic",),
            )
            k == 1 && legend_below!(fig, ax)
        end
        save_figure("graph_matching_accuracy.svg", fig)
    end
end

function plot_edit_distance()
    rows = read_csv("edit_distance.csv")
    with_theme(THEME) do
        fig = Figure()
        ax = seconds_axis(fig[1, 1]; xlabel="number of vertices n", width=600)
        order = [string.(nameof.(GED_FORMULATIONS)); "NetworkX"]
        series!(ax, rows; x=r -> r.n, y=r -> r.time, group=r -> r.method, order)
        ax.xticks = sort(unique(r.n for r in rows))
        legend_below!(fig, ax; nbanks=2)
        save_figure("edit_distance_time.svg", fig)
    end
end

function plot_paths_flows()
    rows = read_csv("paths_flows.csv")
    with_theme(THEME) do
        fig = Figure()
        for (k, problem) in enumerate(("shortest path", "min-cost flow", "assignment"))
            ax = seconds_axis(
                fig[1, k]; title=problem, xscale=log10, xlabel="number of vertices n"
            )
            series!(
                ax,
                filter(r -> r.problem == problem, rows);
                x=r -> r.n,
                y=r -> r.time,
                group=r -> r.library,
            )
            Legend(fig[2, k], ax; framevisible=false, tellwidth=false, valign=:top)
        end
        save_figure("paths_flows_time.svg", fig)
    end
end

function plot_cliques()
    rows = read_csv("cliques.csv")
    problems = ("maximum weight clique", "maximum independent set", "minimum vertex cover")
    with_theme(THEME) do
        fig = Figure()
        for (k, problem) in enumerate(problems)
            ax = seconds_axis(
                fig[1, k]; title=problem, xscale=log10, xlabel="number of vertices n"
            )
            series!(
                ax,
                filter(r -> r.problem == problem, rows);
                x=r -> r.n,
                y=r -> r.time,
                group=r -> r.method,
            )
            Legend(fig[2, k], ax; framevisible=false, tellwidth=false, valign=:top)
        end
        save_figure("cliques_time.svg", fig)

        # Quality of the heuristics, relative to the optimum found by the exact methods
        fig = Figure()
        for (k, problem) in enumerate(problems[2:3])
            ax = Axis(
                fig[1, k];
                title=problem,
                xscale=log10,
                xlabel="number of vertices n",
                ylabel="size / optimal size",
            )
            prows = filter(r -> r.problem == problem && r.valid == "true", rows)
            relative = []
            for n in unique(r.n for r in prows)
                exact = [r.value for r in prows if r.n == n && r.exact == "true"]
                isempty(exact) && continue
                opt = problem == "maximum independent set" ? maximum(exact) : minimum(exact)
                for r in prows
                    r.n == n && push!(relative, (; n, r.method, ratio=r.value / opt))
                end
            end
            series!(
                ax,
                relative;
                x=r -> r.n,
                y=r -> r.ratio,
                group=r -> r.method,
                halo=("NetworkX",),
            )
            Legend(fig[2, k], ax; framevisible=false, tellwidth=false, valign=:top)
        end
        save_figure("cliques_quality.svg", fig)
    end
end
