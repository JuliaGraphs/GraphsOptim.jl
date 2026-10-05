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
    x = @variable(
        model,
        [i in 1:n, j in 1:m];
        lower_bound=0,
        integer=integer,
        base_name=String(var_name)
    )
    model[Symbol(var_name)] = x
    @objective(model, Min, dot(edge_cost, x))
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

"""
    linear_sum_assignment(cost; maximize=false)

Solve the square linear assignment problem with a shortest augmenting path algorithm, without relying on a JuMP optimizer.

Return a vector `σ` such that row `i` is assigned to column `σ[i]`, minimizing (or maximizing) `∑ᵢ cost[i, σ[i]]`.

# References

- Crouse, "On implementing 2D rectangular assignment algorithms", IEEE Transactions on Aerospace and Electronic Systems (2016), also used by `scipy.optimize.linear_sum_assignment`
"""
function linear_sum_assignment(cost::AbstractMatrix{<:Real}; maximize::Bool=false)
    n = size(cost, 1)
    size(cost, 2) == n || throw(DimensionMismatch("the cost matrix must be square"))
    T = float(eltype(cost))
    s = maximize ? -one(T) : one(T)
    # The columns of `cost` play the role of workers and its rows the role of jobs,
    # so that the innermost loop runs down a column (memory-contiguous in Julia).
    u = zeros(T, n)  # dual variables of the workers
    v = zeros(T, n)  # dual variables of the jobs
    shortest = Vector{T}(undef, n)
    path = zeros(Int, n)
    job4worker = zeros(Int, n)
    worker4job = zeros(Int, n)
    visited_workers = Vector{Bool}(undef, n)
    visited_jobs = Vector{Bool}(undef, n)
    remaining = Vector{Int}(undef, n)
    @inbounds for current in 1:n
        fill!(shortest, typemax(T))
        fill!(visited_workers, false)
        fill!(visited_jobs, false)
        for k in 1:n
            remaining[k] = n - k + 1
        end
        num_remaining = n
        min_val = zero(T)
        w = current
        sink = 0
        # Find a shortest augmenting path starting from worker `current`
        while sink == 0
            visited_workers[w] = true
            index = 0
            lowest = typemax(T)
            uw = u[w]
            for k in 1:num_remaining
                j = remaining[k]
                r = min_val + s * cost[j, w] - uw - v[j]
                if r < shortest[j]
                    path[j] = w
                    shortest[j] = r
                end
                if shortest[j] < lowest || (shortest[j] == lowest && worker4job[j] == 0)
                    lowest = shortest[j]
                    index = k
                end
            end
            min_val = lowest
            isfinite(min_val) || throw(ArgumentError("the cost matrix is infeasible"))
            j = remaining[index]
            if worker4job[j] == 0
                sink = j
            else
                w = worker4job[j]
            end
            visited_jobs[j] = true
            remaining[index] = remaining[num_remaining]
            num_remaining -= 1
        end
        # Update the dual variables
        u[current] += min_val
        for w in 1:n
            if visited_workers[w] && w != current
                u[w] += min_val - shortest[job4worker[w]]
            end
        end
        for j in 1:n
            if visited_jobs[j]
                v[j] -= min_val - shortest[j]
            end
        end
        # Augment the assignment along the path
        j = sink
        while true
            w = path[j]
            worker4job[j] = w
            job4worker[w], j = j, job4worker[w]
            w == current && break
        end
    end
    return worker4job
end
