_graphs_not_isomorphic() = throw(ArgumentError("The graphs are not isomorphic"))

"""
    graph_isomorphim(g, h; optimizer=HiGHS.Optimizer)

Return an isomorphism from `g` to `h` as a dictionary mapping each vertex of `g`
to a vertex of `h`. Throw an error if the graphs are not isomorphic.

The returned dictionary `f` preserves adjacency and, for graph types whose
`adjacency_matrix` contains edge weights, preserves those weights as well. Both
directed and undirected graphs are supported, but a directed graph is never
considered isomorphic to an undirected graph.

# Supported graph types

- undirected graphs, such as `SimpleGraph`;
- directed graphs, such as `SimpleDiGraph`;
- weighted graphs, such as `SimpleWeightedGraph` and `SimpleWeightedDiGraph`.

For weighted graphs, an isomorphism must preserve both adjacency and edge
weights.

# Keyword arguments

- `optimizer`: JuMP-compatible mixed-integer solver (default is `HiGHS.Optimizer`)

# Examples

```jldoctest
julia> using Graphs, GraphsOptim

julia> f = graph_isomorphim(path_graph(3), SimpleGraph(Edge.([(2, 3), (3, 1)])));

julia> all(has_edge(path_graph(3), u, v) == has_edge(SimpleGraph(Edge.([(2, 3), (3, 1)])), f[u], f[v]) for u in 1:3, v in 1:3)
true
```
"""
function graph_isomorphim(
    g::AbstractGraph, h::AbstractGraph; optimizer=HiGHS.Optimizer
)::Dict{Int,Int}
    n = nv(g)
    if n != nv(h) || ne(g) != ne(h) || is_directed(g) != is_directed(h)
        _graphs_not_isomorphic()
    end
    n == 0 && return Dict{Int,Int}()

    A = adjacency_matrix(g)
    B = adjacency_matrix(h)

    # These inexpensive invariants avoid building a MIP for many non-isomorphic pairs.
    if !isapprox(sort(vec(sum(A; dims=1))), sort(vec(sum(B; dims=1)))) ||
        !isapprox(sort(vec(sum(A; dims=2))), sort(vec(sum(B; dims=2))))
        _graphs_not_isomorphic()
    end

    model = Model(optimizer)
    set_silent(model)
    @variable(model, X[1:n, 1:n], Bin)

    @constraint(model, [i in 1:n], sum(X[i, j] for j in 1:n) == 1)
    @constraint(model, [j in 1:n], sum(X[i, j] for i in 1:n) == 1)
    @constraint(
        model,
        [i in 1:n, k in 1:n],
        sum(A[i, j] * X[j, k] for j in 1:n) == sum(X[i, j] * B[j, k] for j in 1:n),
    )

    optimize!(model)
    status = termination_status(model)
    if status == INFEASIBLE
        _graphs_not_isomorphic()
    elseif status != OPTIMAL
        error("The graph isomorphism solver terminated with status $status")
    end

    return Dict(i => findfirst(j -> value(X[i, j]) > 0.5, 1:n) for i in 1:n)
end

# Algorithm contributed by Ed Scheinerman.
