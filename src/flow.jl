"""
    min_cost_flow!(
        model,
        g, vertex_demand, edge_cost, edge_min_capacity, edge_max_capacity;
        var_name, integer
    )

Modify a JuMP model by adding the variable, constraints and objective necessary to compute a minimum cost flow over a directed graph.

The flow variable will be named `var_name`, see [`min_cost_flow`](@ref) for details on the other arguments.
"""
function min_cost_flow!(
    model::Model,
    g::AbstractGraph,
    vertex_demand::AbstractVector,
    edge_cost::AbstractMatrix,
    edge_min_capacity::AbstractMatrix,
    edge_max_capacity::AbstractMatrix;
    var_name,
    integer::Bool,
)
    if is_directed(g)
        edge_tuples = [(src(ed), dst(ed)) for ed in edges(g)]
    else
        edge_tuples = vcat(
            [(src(ed), dst(ed)) for ed in edges(g)], [(dst(ed), src(ed)) for ed in edges(g)]
        )
    end

    f = @variable(model, [edge_tuples]; integer=integer, base_name=String(var_name))
    model[Symbol(var_name)] = f

    for (u, v) in edge_tuples
        if isfinite(edge_min_capacity[u, v])
            @constraint(model, f[(u, v)] >= edge_min_capacity[u, v])
        end
        if isfinite(edge_max_capacity[u, v])
            @constraint(model, f[(u, v)] <= edge_max_capacity[u, v])
        end
    end

    for v in vertices(g)
        inflow, outflow = AffExpr(0), AffExpr(0)
        for u in inneighbors(g, v)
            inflow += f[(u, v)]
        end
        for w in outneighbors(g, v)
            outflow += f[(v, w)]
        end
        @constraint(model, inflow == vertex_demand[v] + outflow)
    end

    obj = objective_function(model)
    for (u, v) in edge_tuples
        add_to_expression!(obj, edge_cost[u, v], f[(u, v)])
    end
    @objective(model, Min, obj)

    return model
end

"""
    min_cost_flow(
        g, vertex_demand, edge_cost, edge_min_capacity, edge_max_capacity;
        integer, optimizer
    )

Compute a minimum cost flow over a directed graph.

# Arguments

- `g::Graphs.AbstractGraph`: a directed graph `G = (V, E)`
- `vertex_demand::AbstractVector`: a vector in `Rⱽ` giving the flow requested by each vertex (should be positive for sinks, negative for sources and zero elsewhere)
- `edge_cost::AbstractMatrix`: a vector in `Rᴱ` giving the cost of a unit of flow on each edge
- `edge_min_capacity::AbstractMatrix`: a vector in `Rᴱ` giving the minimum flow allowed on each edge
- `edge_max_capacity::AbstractMatrix`: a vector in `Rᴱ` giving the maximum flow allowed on each edge

# Keyword arguments

- `integer::Bool`: whether the flow should be integer-valued or real-valued
- `optimizer`: JuMP-compatible solver (default is `HiGHS.Optimizer`)
"""
function min_cost_flow(
    g,
    vertex_demand,
    edge_cost,
    edge_min_capacity=Zeros(nv(g), nv(g)),
    edge_max_capacity=Fill(Inf, nv(g), nv(g));
    integer=false,
    optimizer=HiGHS.Optimizer,
)
    model = Model(optimizer)
    set_silent(model)
    min_cost_flow!(
        model,
        g,
        vertex_demand,
        edge_cost,
        edge_min_capacity,
        edge_max_capacity;
        var_name=:flow,
        integer,
    )
    optimize!(model)
    @assert termination_status(model) == OPTIMAL

    flow_values = value.(model[:flow])
    I = src.(edges(g))
    J = dst.(edges(g))
    V = [flow_values[(i, j)] for (i, j) in zip(I, J)]
    flow_matrix = sparse(I, J, V, nv(g), nv(g))
    return flow_matrix
end
