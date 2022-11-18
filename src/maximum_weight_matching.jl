"""
    MatchingResult{U}

A type representing the result of a matching algorithm.

# Fields

- `weight::U`: total weight of the matching
- `mate::Vector{Int}`: pairwise assignment.

`mate[i] = j` if vertex `i` is matched to vertex `j`, and `mate[i] = -1` for unmatched vertices.
"""
struct MatchingResult{U<:Real}
    weight::U
    mate::Vector{Int}
end

## Maximum weight matching

"""
    maximum_weight_matching(g, w; optimizer)

Given a graph `g` and an matrix `w` of edge weights, return a matching  ([`MatchingResult`](@ref) object) with the maximum total weight.

If no weight matrix is given, all edges will be considered to have weight 1
(results in max cardinality matching).
Edges in `g` that are not present in `w` will not be considered for the matching.

A JuMP-compatible solver must be provided with the `optimizer` argument.

The efficiency of the algorithm depends on the input graph:
  - If the graph is bipartite, then the LP relaxation is integral.
  - If the graph is not bipartite, then it requires a MIP solver and the computation time may grow exponentially.
"""
function maximum_weight_matching(
    g::Graph, w::AbstractMatrix{U}=default_weights(g); optimizer
) where {U<:Real}
    model = Model(optimizer)
    n = nv(g)
    edge_list = collect(edges(g))

    # put the edge weights in w in the right order to be compatible with edge_list
    for j in 1:n
        for i in 1:n
            if i > j && w[i, j] > zero(U) && w[j, i] < w[i, j]
                w[j, i] = w[i, j]
            end
            if Edge(i, j) ∉ edge_list
                w[i, j] = zero(U)
            end
        end
    end

    if is_bipartite(g)
        @variable(model, x[edge_list] >= 0) # no need to enforce integrality
    else
        @variable(model, x[edge_list] >= 0, Int) # requires MIP solver
    end
    @objective(model, Max, sum(x[e] * w[src(e), dst(e)] for e in edge_list))

    @constraint(model, c1[i=1:n], sum(x[Edge(minmax(i, j))] for j in neighbors(g, i)) <= 1)
    optimize!(model)
    status = JuMP.termination_status(model)
    status != MOI.OPTIMAL && error("JuMP solver failed to find optimal solution.")
    solution = value.(x)
    cost = objective_value(model)

    mate = fill(-1, n)
    for e in edge_list
        if solution[e] >= 1 - 1e-5 # Some tolerance to numerical approximations by the solver.
            mate[src(e)] = dst(e)
            mate[dst(e)] = src(e)
        end
    end

    return MatchingResult(cost, mate)
end
