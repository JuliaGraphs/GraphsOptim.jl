# Graph edit distance: the ILP formulations of GraphsOptim.jl against the exact
# algorithm of NetworkX, with uniform edit costs.

const GED_FORMULATIONS = [F1, F1prime, F1plus, F2, F2minus, F2plus, FORI]

"""
    edited_graph(rng, g, k)

Return a copy of `g` where `k` random vertex pairs have been flipped (edge added or removed), and the vertices shuffled.
"""
function edited_graph(rng::AbstractRNG, g::SimpleGraph, k::Integer)
    h = copy(g)
    pairs = [(i, j) for j in vertices(g) for i in 1:(j - 1)]
    for (i, j) in pairs[randperm(rng, length(pairs))[1:k]]
        has_edge(h, i, j) ? rem_edge!(h, i, j) : add_edge!(h, i, j)
    end
    return h[randperm(rng, nv(h))]
end

function run_edit_distance()
    nx = pyimport("networkx")
    # compilation
    let g = erdos_renyi(4, 0.5; rng=Xoshiro(0))
        for formulation in GED_FORMULATIONS
            GraphsOptim.edit_distance(
                g, edited_graph(Xoshiro(0), g, 1); formulation, optimizer=OPTIMIZER
            )
        end
        nx.graph_edit_distance(to_networkx(g), to_networkx(edited_graph(Xoshiro(0), g, 1)))
    end

    rows = []
    timed_out = Set{String}()
    for n in 4:2:24
        rng = Xoshiro(n)
        G = erdos_renyi(n, 0.4; rng)
        H = edited_graph(rng, G, max(2, n ÷ 3))
        for formulation in GED_FORMULATIONS
            method = string(nameof(formulation))
            method in timed_out && continue
            t, res = measure(
                () -> GraphsOptim.edit_distance(G, H; formulation, optimizer=OPTIMIZER);
                min_time=0,
            )
            too_slow(t) && push!(timed_out, method)
            value = isnothing(res) ? NaN : res.objective_value
            push!(rows, (; n, edges_G=ne(G), edges_H=ne(H), method, time=t, value))
        end
        if !("NetworkX" in timed_out)
            t, value = pymeasure(
                nx.graph_edit_distance,
                to_networkx(G),
                to_networkx(H);
                timeout=TIME_LIMIT,
                min_time=0,
            )
            # On timeout, NetworkX returns the best value found so far, which is not proven optimal
            value = t >= TIME_LIMIT ? NaN : pyconvert(Float64, value)
            t = t >= TIME_LIMIT ? NaN : t
            too_slow(t) && push!(timed_out, "NetworkX")
            push!(
                rows, (; n, edges_G=ne(G), edges_H=ne(H), method="NetworkX", time=t, value)
            )
        end
        @info "Edit distance" n
    end
    write_csv("edit_distance.csv", rows)
    return nothing
end
