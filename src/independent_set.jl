
"""
    maximum_independent_set!(model, g; var_name)

Computes in-place in the JuMP model a maximum-weighted independent set of `g`.
An optional `vertex_weights` vector can be passed to the graph, defaulting to uniform weights (computing a maximum size independent set).
"""
function maximum_weight_independent_set!(
    model::Model, g::AbstractGraph; binary::Bool=true, var_name, vertex_weights=ones(nv(g))
)
    if is_directed(g)
        throw(ArgumentError("The graph must not be directed"))
    end
    g_vertices = collect(vertices(g))
    f = @variable(model, [g_vertices]; binary=binary, base_name=String(var_name))
    model[Symbol(var_name)] = f
    @constraint(
        model,
        covering_constraint[i=1:nv(g), j=1:nv(g); i ≠ j && has_edge(g, i, j)],
        f[i] + f[j] <= 1,
    )
    obj = objective_function(model)
    add_to_expression!(obj, dot(f, vertex_weights))
    @objective(model, Max, obj)
    return model
end

"""
    maximum_weight_independent_set(g; optimizer, binary, vertex_weights)

Computes a maximum-weighted independent set or stable set of `g`.
"""
function maximum_weight_independent_set(
    g::AbstractGraph;
    binary::Bool=true,
    vertex_weights=ones(nv(g)),
    optimizer=HiGHS.Optimizer,
)
    model = Model(optimizer)
    set_silent(model)
    maximum_weight_independent_set!(
        model, g; binary=binary, vertex_weights=vertex_weights, var_name=:stable
    )
    optimize!(model)
    @assert termination_status(model) == OPTIMAL
    stable_variables = Vector(model[:stable])
    stable_vertices = findall(v -> value(v) > 0.5, stable_variables)
    return stable_vertices
end
