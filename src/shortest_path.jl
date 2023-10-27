"""
	shortest_path!(
		model,
		g, source, target, edge_cost;
        var_name, integer
	)

Modify a JuMP model by adding the variable, constraints and objective necessary to compute the shortest path between 2 vertices over a graph.

The edge selection variable will be named `var_name`, see [`shortest_path`](@ref) for details on the other arguments.
"""
function shortest_path!(
    model::Model,
    g::AbstractGraph,
    source::Int,
    target::Int,
    edge_cost::AbstractMatrix;
    var_name,
    integer::Bool=true,
)
    vertex_demand = zeros(Int, nv(g))
    vertex_demand[source] = -1
    vertex_demand[target] = +1
    edge_min_capacity = Zeros(nv(g), nv(g))
    edge_max_capacity = Ones(nv(g), nv(g))
    min_cost_flow!(
        model,
        g,
        vertex_demand,
        edge_cost,
        edge_min_capacity,
        edge_max_capacity;
        var_name,
        integer,
    )
    return model
end

"""
	shortest_path(
		g, source, target, edge_cost;
        integer, optimizer
	)

Compute the shortest path from a source to a target vertex over a graph.
Returns a sequence of vertices from source to destination.

# Arguments

- `g::Graphs.AbstractGraph`: a directed or undirected graph `G = (V, E)`
- `source::Int`: source vertex of the path
- `target::Int`: target vertex of the path
- `edge_cost::AbstractMatrix`: a vector in `Rᴱ` giving the cost of traversing the edge

# Keyword arguments

- `integer::Bool`: whether the path should be integer-valued or real-valued (default is `true`)
- `optimizer`: JuMP-compatible solver (default is `HiGHS.Optimizer`)
"""
function shortest_path(
    g::AbstractGraph,
    source::Int,
    target::Int,
    edge_cost::AbstractMatrix;
    integer::Bool=true,
    optimizer=HiGHS.Optimizer,
)
    model = Model(optimizer)
    set_silent(model)

    shortest_path!(model, g, source, target, edge_cost; var_name=:path, integer=integer)
    optimize!(model)
    @assert termination_status(model) == OPTIMAL

    edge_selection = value.(model[:path])
    valid_edges = [p for p in first(axes(edge_selection)) if edge_selection[p] > 0]

    # Parse valid path indices to actual path (sequence of visited vertices)
    path = ones(typeof(source), length(valid_edges) + 1)
    path[begin] = source

    for step in 1:length(valid_edges)
        prev = path[step]
        path[step + 1] = first(edge[end] for edge in valid_edges if edge[begin] == prev)
    end

    return path
end
