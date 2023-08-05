function indvec(s, n)
    x = zeros(n)
    x[s] .= 1
    return x
end

"""
    fractional_chromatic_number(g; optimizer)

Compute the fractional chromatic number of a graph.  Gives the same result as
fractional_clique_number, though one function may run faster than the other.
Beware: this can run very slowly for graphs of any substantial size.

# Keyword arguments

- `optimizer`: JuMP-compatible solver (default is `GLPK.Optimizer`)
"""
function fractional_chromatic_number(
    g::AbstractGraph{T}, optimizer=GLPK.Optimizer
) where {T<:Integer}
    ss = maximal_cliques(complement(g))
    M = hcat(indvec.(ss, nv(g))...)

    model = Model(optimizer)
    set_silent(model)
    @variable(model, x[1:length(ss)] >= 0)
    @constraint(model, M * x .>= 1)
    @objective(model, Min, sum(x))
    optimize!(model)
    return objective_value(model)
end

"""
    fractional_clique_number(g; optimizer)

Compute the fractional clique number of a graph.  Gives the same result as
fractional_chromatic_number, though one function may run faster than the other.
Beware: this can run very slowly for graphs of any substantial size.

# Keyword arguments

- `optimizer`: JuMP-compatible solver (default is `GLPK.Optimizer`)
"""
function fractional_clique_number(
    g::AbstractGraph{T}, optimizer=GLPK.Optimizer
) where {T<:Integer}
    model = Model(optimizer)
    set_silent(model)
    @variable(model, x[1:nv(g)] >= 0)
    for clique in maximal_cliques(complement(g))
        @constraint(model, sum(x[clique]) <= 1)
    end
    @objective(model, Max, sum(x))
    optimize!(model)
    return objective_value(model)
end
