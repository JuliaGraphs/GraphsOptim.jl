"""
    minimum_cost_flow!(
        model, g, node_demand, edge_cost, edge_min_capacity, edge_max_capacity; flowvar_name
    )

Compute a flow over a directed graph `g` that satisfies the demand at each node and capacity constraints along each edge, while minimizing the total edge cost.

The first argument is a JuMP `model` that will be modified with the addition of: 1) a new variable called `flowvar_name`, 2) flow constraints and 3) an additional objective term to minimize.
"""
@traitfn function minimum_cost_flow!(
    model::Model,
    g::G,
    node_demand::AbstractVector,
    edge_cost::AbstractMatrix,
    edge_min_capacity::AbstractMatrix=Zeros(nv(g), nv(g)),
    edge_max_capacity::AbstractMatrix=Fill(Inf, nv(g), nv(g));
    flowvar_name::Symbol=:flow,
) where {G <: AbstractGraph; IsDirected{G}}
    edge_tuples = [(src(ed), dst(ed)) for ed in edges(g)]

    f = @variable(model, [edge_tuples]; base_name=String(flowvar_name))
    model[flowvar_name] = f

    for (u, v) in edge_tuples
        @constraint(model, f[(u, v)] >= edge_min_capacity[u, v])
        if isfinite(edge_max_capacity[u, v])
            @constraint(model, f[(u, v)] <= edge_max_capacity[u, v])
        end
    end

    for v in vertices(g)
        inflow = sum(f[(u, v)] for u in inneighbors(g, v))
        outflow = sum(f[(v, w)] for w in outneighbors(g, v))
        @constraint(model, inflow == node_demand[v] + outflow)
    end

    obj = objective_function(model)
    for (u, v) in edge_tuples
        add_to_expression!(obj, edge_cost[u, v], f[(u, v)])
    end
    @objective(model, Min, obj)

    return model
end
