"""
	shortest_path!(
		model,
		g, source, target, edge_cost, var_name
	)

Modify a JuMP model by adding the variable, constraints and objective necessary to compute the shortest path between 2 vertices over a graph

The edge selection variable will be named `var_name`
"""
function shortest_path!(
    model::Model,
    g::AbstractGraph,
    source::Int,
    target::Int,
    edge_cost::AbstractMatrix;
    var_name,
)
    dg = g
    if !is_directed(g)
        dg = DiGraph(g)
    end
    edge_tuples = [(src(ed), dst(ed)) for ed in edges(dg)]

    # 1 if an edge if traversed and 0 otherwise
    traversed = @variable(model, 0 <= x[e=edge_tuples] <= 1; integer=true)
    model[Symbol(var_name)] = traversed

    # Exit from source
    @constraint(
        model,
        sum(traversed[(source, next_v)] for next_v in outneighbors(dg, source)) -
        sum(traversed[(prev_v, source)] for prev_v in inneighbors(dg, source)) == 1
    )
    # Arrive at destination
    @constraint(
        model,
        sum(traversed[(target, next_v)] for next_v in outneighbors(dg, target)) +
        sum(traversed[(prev_v, target)] for prev_v in inneighbors(dg, target)) == 1
    )
    # Go out from a vertex is only possible by going in first
    for u in vertices(dg)
        if u in (source, target)
            continue
        end
        @constraint(
            model,
            sum(traversed[(u, next_v)] for next_v in outneighbors(dg, u)) -
            sum(traversed[(prev_v, u)] for prev_v in inneighbors(dg, u)) == 0
        )
    end

    @objective(
        model, Min, sum(edge_cost[u, v] * traversed[(u, v)] for (u, v) in edge_tuples)
    )

    return model
end

"""
	shortest_path(
		g, source, target, edge_cost, optimizer
	)

Compute the shortest path from a source to a target vertex over a graph.
Returns a sequence of vertices from source to destination.

# Arguments

- `g::Graphs.AbstractGraph`: a directed or undirected graph `G = (V, E)`
- `source::Int`: source vertex of the path
- `target::Int`: target/destination vertex of the path
- `edge_cost::AbstractMatrix`: a vector in `Rᴱ` giving the cost of traversing the edge

# Keyword arguments

- `optimizer`: JuMP-compatible solver (default is `HiGHS.Optimizer`)
"""
function shortest_path(
    g::AbstractGraph,
    source::Int,
    target::Int,
    edge_cost::AbstractMatrix;
    optimizer=HiGHS.Optimizer,
)
    model = Model(optimizer)
    set_silent(model)

    shortest_path!(model, g, source, target, edge_cost; var_name=:path)
    optimize!(model)
    @assert termination_status(model) == OPTIMAL

    edge_selection = value.(model[:path])
    valid_edges = [p for p in first(axes(edge_selection)) if edge_selection[p] > 0]

    # Parse valid path indices to actual path (sequence of traversed vertices)
    path = ones(typeof(source), length(valid_edges) + 1)
    path[begin] = source

    for step in 1:length(valid_edges)
        prev = path[step]
        edge_index = findfirst(edge[begin] == prev for edge in valid_edges)
        path[step + 1] = valid_edges[edge_index][end]
    end

    return path
end