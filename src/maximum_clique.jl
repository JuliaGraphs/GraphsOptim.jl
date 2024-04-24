
"""
    maximum_weight_clique!(model, g; var_name)

Computes in-place in the JuMP model a maximum-weighted clique of `g`.
An optional `vertex_weights` vector can be passed to the graph, defaulting to uniform weights (computing a maximum size clique).
"""
function maximum_weight_clique!(
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
        packing_constraint[i=1:nv(g), j=1:nv(g); i ≠ j && !has_edge(g, i, j)],
        f[i] + f[j] <= 1,
    )
    @objective(model, Max, dot(f, vertex_weights))
    return model
end

"""
    maximum_weight_clique(g; optimizer, binary, vertex_weights)

Computes a maximum-weighted clique of `g`.
"""
function maximum_weight_clique(
    g::AbstractGraph;
    binary::Bool=true,
    vertex_weights=ones(nv(g)),
    optimizer=HiGHS.Optimizer,
)
    model = Model(optimizer)
    set_silent(model)
    maximum_weight_clique!(
        model, g; binary=binary, vertex_weights=vertex_weights, var_name=:clique
    )
    optimize!(model)
    @assert termination_status(model) == OPTIMAL
    clique_variables = Vector(model[:clique])
    clique_vertices = findall(v -> value(v) > 0.5, clique_variables)
    return clique_vertices
end
