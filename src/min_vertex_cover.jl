"""
    min_vertex_cover!(
        model, g;
        var_name
    )

Modify a JuMP model by adding the variable, constraints and objective necessary to compute a minimum vertex cove of an undirected graph.

The vertex cover indicator variable will be named `var_name`
"""
function min_vertex_cover!(model::Model, g::AbstractGraph; integer::Bool, var_name)
    if is_directed(g)
        throw(ArgumentError("The graph must not be directed"))
    end

    ver = collect(vertices(g))

    f = @variable(model, [ver]; integer=integer, base_name=String(var_name))
    model[Symbol(var_name)] = f

    for v in ver
        @constraint(model, f[v] >= 0)
        @constraint(model, f[v] <= 1)
    end

    for e in edges(g)
        @constraint(model, f[src(e)] + f[dst(e)] >= 1)
    end

    obj = objective_function(model)
    for v in ver
        add_to_expression!(obj, f[v])
    end
    @objective(model, Min, obj)

    return model
end

"""
    min_vertex_cover(
        g, optimizer
    )

Compute a minimum vertex cover of an undirected graph.

# Arguments

- `g::Graphs.AbstractGraph`: an undirected graph `G = (V, E)`

# Keyword arguments

- `optimizer`: JuMP-compatible solver (default is `HiGHS.Optimizer`)
"""
function min_vertex_cover(g::AbstractGraph; integer::Bool=true, optimizer=HiGHS.Optimizer)
    model = Model(optimizer)
    set_silent(model)

    min_vertex_cover!(model, g; integer=integer, var_name=:cover)

    optimize!(model)
    @assert termination_status(model) == OPTIMAL

    cover_values = Vector(value.(model[:cover]))
    cover_vertices = findall(==(1), cover_values)

    return cover_vertices
end
