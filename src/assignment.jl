"""
    min_cost_assignment!(
        model,
        edge_cost;
        var_name, integer
    )

Modify a JuMP model by adding the variable, constraints and objective necessary to compute a minimum cost assignment over a bipartite graph.

The assignment variable will be named `var_name`, see [`min_cost_assignment`](@ref) for details on the other arguments.
"""
function min_cost_assignment!(
    model::Model, edge_cost::AbstractMatrix; var_name, integer::Bool
)
    n, m = size(edge_cost)
    x = @variable(model, [i in 1:n, j in 1:m]; integer=integer, base_name=String(var_name))
    model[Symbol(var_name)] = x
    @objective(model, Min, dot(edge_cost, x))
    @constraint(model, x .>= 0)
    @constraint(model, [i in 1:n], sum(x[i, :]) == 1)
    @constraint(model, [j in 1:m], sum(x[:, j]) == 1)
    return nothing
end

"""
    min_cost_assignment(
        edge_cost;
        optimizer
    )

Compute a minimum cost assignment over a bipartite graph.

# Arguments

- `edge_cost::AbstractMatrix`: a matrix in `Rᵁˣⱽ` giving the cost of matching `u ∈ U` to `v ∈ V`

# Keyword arguments

- `integer::Bool`: whether the flow should be integer-valued or real-valued
- `optimizer`: JuMP-compatible solver (default is `HiGHS.Optimizer`)
"""
function min_cost_assignment(
    edge_cost::AbstractMatrix; integer::Bool, optimizer=HiGHS.Optimizer
)
    model = Model(optimizer)
    set_silent(model)
    min_cost_assignment!(model, edge_cost; var_name=:assignment, integer=integer)
    optimize!(model)
    @assert termination_status(model) == OPTIMAL
    return value.(model[:assignment])
end
