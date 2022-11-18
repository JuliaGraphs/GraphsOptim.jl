
## Maximum weight maximal matching

"""
    maximum_weight_maximal_matching(g, w[, cutoff]; optimizer)

Given a bipartite graph `g` and a matrix `w` of edge weights, return a matching ([`MatchingResult`](@ref) object) with the maximum total weight among the ones that contain the largest number of edges.

If no weight matrix is given, all edges will be considered to have weight 1
(results in max cardinality matching).
Edges in `g` that are not present in `w` will not be considered for the matching.

A JuMP-compatible solver must be provided with the `optimizer` argument.

A `cutoff` argument can be given to reduce computation time by excluding edges with weights lower than the specified value.

The algorithm relies on a linear relaxation on of the matching problem, which is
guaranteed to have integer solution on bipartite graphs.
"""
function maximum_weight_maximal_matching(
    g::Graph, w::AbstractMatrix{T}=default_weights(g); optimizer
) where {T<:Real}
    # TODO support for graphs with zero degree nodes
    # TODO apply separately on each connected component
    bpmap = bipartite_map(g)
    length(bpmap) != nv(g) && error("Graph is not bipartite")
    v1 = findall(isequal(1), bpmap)
    v2 = findall(isequal(2), bpmap)
    if length(v1) > length(v2)
        v1, v2 = v2, v1
    end

    nedg = 0
    edgemap = Dict{Edge,Int}()

    for j in 1:size(w, 2)
        for i in 1:size(w, 1)
            if w[i, j] > 0.0
                nedg += 1
                edgemap[Edge(i, j)] = nedg
                edgemap[Edge(j, i)] = nedg
            end
        end
    end

    model = Model(optimizer)
    @variable(model, x[1:length(w)] >= 0)

    for i in v1
        idx = Vector{Int}()
        for j in neighbors(g, i)
            if haskey(edgemap, Edge(i, j))
                push!(idx, edgemap[Edge(i, j)])
            end
        end
        if length(idx) > 0
            @constraint(model, sum(x[id] for id in idx) == 1)
        end
    end

    for j in v2
        idx = Vector{Int}()
        for i in neighbors(g, j)
            if haskey(edgemap, Edge(i, j))
                push!(idx, edgemap[Edge(i, j)])
            end
        end

        if length(idx) > 0
            @constraint(model, sum(x[id] for id in idx) <= 1)
        end
    end

    @objective(model, Max, sum(w[src(e), dst(e)] * x[edgemap[e]] for e in keys(edgemap)))

    optimize!(model)
    status = JuMP.termination_status(model)
    status != MOI.OPTIMAL && error("JuMP solver failed to find optimal solution.")
    sol = JuMP.value.(x)

    all(Bool[s == 1 || s == 0 for s in sol]) || error("Found non-integer solution.")

    cost = JuMP.objective_value(model)

    mate = fill(-1, nv(g))
    for e in edges(g)
        if w[src(e), dst(e)] > zero(T)
            inmatch = convert(Bool, sol[edgemap[e]])
            if inmatch
                mate[src(e)] = dst(e)
                mate[dst(e)] = src(e)
            end
        end
    end

    return MatchingResult(cost, mate)
end

"""
    cutoff_weights(w, cutoff)

Copy the weights matrix `w` with all elements below `cutoff` set to 0.
"""
function cutoff_weights(w::AbstractMatrix{T}, cutoff::R) where {T<:Real,R<:Real}
    wnew = copy(w)
    for j in 1:size(w, 2)
        for i in 1:size(w, 1)
            if wnew[i, j] < cutoff
                wnew[i, j] = zero(T)
            end
        end
    end
    return wnew
end

function maximum_weight_maximal_matching(
    g::Graph, w::AbstractMatrix{T}, cutoff::R; optimizer
) where {T<:Real,R<:Real}
    return maximum_weight_maximal_matching(g, cutoff_weights(w, cutoff); optimizer)
end
