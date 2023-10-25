function vertex_cover!(
    model::Model,
    g::AbstractGraph;
    var_name,
)
    if is_directed(g)
        throw(ArgumentError("The graph must not be directed"))
    end
    
    ver = collect(vertices(g))
    edge_tuples = [(src(ed), dst(ed)) for ed in edges(g)]

    f = @variable(model, [ver]; integer=true, base_name=String(var_name))
    model[Symbol(var_name)] = f

    for v in ver
        @constraint(model, f[v] >= 0)
        @constraint(model, f[v] <= 1)
    end

    for (u,v) in edge_tuples
        @constraint(model, f[u] + f[v] >= 1)
    end

    obj = objective_function(model)
    for v in ver
        add_to_expression!(obj, f[v])
    end
    @objective(model, Min, obj)

    return model
end



function vertex_cover(
    g::AbstractGraph,
    optimizer=HiGHS.Optimizer,
)
    model = Model(optimizer)
    set_silent(model)

    vertex_cover!(
        model,
        g;
        var_name=:cover,
    )

    optimize!(model)
    @assert termination_status(model) == OPTIMAL

    cover_values = Vector(value.(model[:cover]))
    cover_vertices = findall(==(1), cover_values)

    return cover_vertices
end
