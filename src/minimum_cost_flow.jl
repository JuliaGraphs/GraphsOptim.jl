@doc raw"""
    minimum_cost_flow!(
        model, g, vertex_demand, edge_cost, edge_min_capacity, edge_max_capacity;
        flowvar_name, integer
    )

Modify a JuMP model by adding the variables, constraints and objective necessary to compute a minimum cost flow over a directed graph.

# Arguments

- `model::JuMP.Model`: the optimization model to modify
- `g::Graphs.AbstractGraph`: a directed graph $G = (V, E)$
- `vertex_demand::AbstractVector`: a vector $d \in \mathbb{R}^V$ giving the flow requested by each vertex (should be positive for sinks, negative for sources and zero elsewhere)
- `edge_cost::AbstractMatrix`: a vector $c \in \mathbb{R}^E$ giving the cost of a unit of flow on each edge
- `edge_min_capacity::AbstractMatrix`: a vector $a \in \mathbb{R}^E$ giving the minimum flow allowed on each edge
- `edge_max_capacity::AbstractMatrix`: a vector $b \in \mathbb{R}^E$ giving the maximum flow allowed on each edge

# Keyword arguments

- `flowvar_name::Symbol`: the name of the optimization variable containing the flow
- `integer::Bool`: whether the flow should be integer-valued or real-valued

# Mathematical formulation

The objective function is
```math
\min_{f \in \mathbb{R}^E} \sum_{(u, v) \in E} c(u, v) f(u, v)
```
The edge capacity constraint dictates that for all $(u, v) \in E$,
```math
a(u, v) \leq f(u, v) \leq b(u, v)
```
The flow conservation constraint with node demand dictates that for all $v \in V$,
```math
f^-(v) = d(v) + f^+(v)
```
where the incoming flow $f^-(v)$ and outgoing flow $f^+(v)$ are defined as
```math
f^-(v) = \sum_{u \in N^-(v)} f(u, v) \quad \text{and} \quad f^+(v) = \sum_{w \in N^+(v)} f(v, w)
```
"""
@traitfn function minimum_cost_flow!(
    model::Model,
    g::G,
    vertex_demand::AbstractVector,
    edge_cost::AbstractMatrix,
    edge_min_capacity::AbstractMatrix,
    edge_max_capacity::AbstractMatrix;
    flowvar_name,
    integer,
) where {G <: AbstractGraph; IsDirected{G}}
    edge_tuples = [(src(ed), dst(ed)) for ed in edges(g)]

    f = @variable(model, [edge_tuples]; integer=integer, base_name=String(flowvar_name))
    model[Symbol(flowvar_name)] = f

    for (u, v) in edge_tuples
        @constraint(model, f[(u, v)] >= edge_min_capacity[u, v])
        if isfinite(edge_max_capacity[u, v])
            @constraint(model, f[(u, v)] <= edge_max_capacity[u, v])
        end
    end

    for v in vertices(g)
        inflow = sum(f[(u, v)] for u in inneighbors(g, v))
        outflow = sum(f[(v, w)] for w in outneighbors(g, v))
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
    minimum_cost_flow(
        g, vertex_demand, edge_cost, edge_min_capacity, edge_max_capacity;
        flowvar_name, integer
    )

Create a JuMP model, call [`minimum_cost_flow!`](@ref) to populate it, optimize it and return the optimal flow as a sparse matrix.
"""
function minimum_cost_flow(
    g,
    vertex_demand,
    edge_cost,
    edge_min_capacity=Zeros(nv(g), nv(g)),
    edge_max_capacity=Fill(Inf, nv(g), nv(g));
    flowvar_name=:flow,
    integer=false,
)
    model = Model(Optimizer)
    set_silent(model)
    minimum_cost_flow!(
        model,
        g,
        vertex_demand,
        edge_cost,
        edge_min_capacity,
        edge_max_capacity;
        flowvar_name,
        integer,
    )
    optimize!(model)
    @assert termination_status(model) == MOI.OPTIMAL

    flow_values = value.(model[Symbol(flowvar_name)])
    I = src.(edges(g))
    J = dst.(edges(g))
    V = [flow_values[(i, j)] for (i, j) in zip(I, J)]
    flow_matrix = sparse(I, J, V, nv(g), nv(g))
    return flow_matrix
end
