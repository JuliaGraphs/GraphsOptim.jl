"""
Control flow types to communicate the ILP formulation to the edit distance code. The
subtypes are the different formulations, with `FORI` being the best performing one. 
"""
abstract type Formulation end
"""
Control flow type signifying the F1 ILP formulation. See also [`Formulation`](@ref).
"""
struct F1 <: Formulation end
"""
Control flow type signifying the F1' ILP formulation. See also [`Formulation`](@ref).
"""
struct F1prime <: Formulation end
"""
Control flow type signifying the F1+ ILP formulation. See also [`Formulation`](@ref).
"""
struct F1plus <: Formulation end
"""
Control flow type signifying the F2- ILP formulation. See also [`Formulation`](@ref).
"""
struct F2minus <: Formulation end
"""
Control flow type signifying the F2 ILP formulation. See also [`Formulation`](@ref).
"""
struct F2 <: Formulation end
"""
Control flow type signifying the F2+ ILP formulation. See also [`Formulation`](@ref).
"""
struct F2plus <: Formulation end
"""
Control flow type signifying the FORI ILP formulation. See also [`Formulation`](@ref).
"""
struct FORI <: Formulation end

"""
Struct maintaining the variable references for formulations [`F2minus`](@ref), [`F2`](@ref)
and [`F2plus`](@ref). Compared to [`FullVariables`](@ref), the deletion variables are
inferred.
"""
struct ReducedVariables
    x::Matrix{VariableRef}
    y
end

"""
Struct maintaining the variable references for formulations [`F1`](@ref), [`F1prime`](@ref)
and [`F1plus`](@ref).
"""
struct FullVariables
    x::Matrix{VariableRef}
    y
    nodeDelG::Vector{VariableRef}
    nodeDelH::Vector{VariableRef}
    edgeDelG
    edgeDelH
end

"""
Struct maintaining the variable references for formulation [`FORI`](@ref). This formulation
implicitly adds two oriented edges for each edge the second graph to improve topology
constraints, which is reflected in the edge map variables.
"""
struct OrientedVariables
    x::Matrix{VariableRef}
    z
end

"""
Type used for functions applicable to both F1 and F2 derived formulations.
"""
Variables = Union{ReducedVariables,FullVariables}

"""
    create_model_vars_reduced!(model, G, H, bidirectional = false)

Adds necessary variables for edit distance computation between `G` and `H` to `model`. If
`bidirectional` is false, for formulations [`F2minus`](@ref), [`F2`](@ref),
[`F2plus`](@ref). If `bidirectional` is true, `H` is bidirected and the variables are used
for the [`FORI`](@ref) formulation.

The variables model a full mapping between nodes and edges respectively.

# Returns
- [`OrientedVariables`](@ref) or [`ReducedVariables`](@ref) containing the variables.
"""
function create_model_vars_reduced!(
    model::Model, G::AbstractGraph, H::AbstractGraph, bidirectional::Bool=false
)
    @variable(model, x[1:nv(G), 1:nv(H)], Bin)

    # Use a sparse indexed edge variable set, to later do partial sums y[i, :, k, l].
    # Edge variables are always oriented to avoid ambiguity.
    @variable(
        model,
        y[
            i in vertices(G),
            j in vertices(G),
            k in vertices(H),
            l in vertices(H);
            has_edge(G, i, j) && has_edge(H, k, l) && i < j && (k < l || bidirectional),
        ],
        Bin
    )

    if bidirectional
        return OrientedVariables(x, y)
    else
        return ReducedVariables(x, y)
    end
end

"""
Convenience function bind. See [`create_model_vars_reduced!`](@ref).
"""
function create_model_vars_bidirectional!(model::Model, G::AbstractGraph, H::AbstractGraph)
    return create_model_vars_reduced!(model, G, H, true)
end

"""
    create_model_vars_full!(model, G, H)

Adds necessary variables for edit distance computation between `G` and `H` to `model` for
formulations [`F1`](@ref), [`F1prime`](@ref) and [`F1plus`](@ref).

In addition to the normal node and edge map variables (see
[`create_model_vars_reduced!`](@ref)), these formulations explicitly have variables to model
nodes and edges being deleted or added.
"""
function create_model_vars_full!(model::Model, G::AbstractGraph, H::AbstractGraph)
    vars = create_model_vars_reduced!(model, G, H)
    @variable(model, nodeDelG[1:nv(G)], Bin)
    @variable(model, nodeDelH[1:nv(H)], Bin)

    @variable(model, edgeDelG[i in vertices(G), j in vertices(G); has_edge(G, i, j)], Bin)
    @variable(model, edgeDelH[k in vertices(H), l in vertices(H); has_edge(H, k, l)], Bin)

    return FullVariables(vars.x, vars.y, nodeDelG, nodeDelH, edgeDelG, edgeDelH)
end

"""
Models arbitrary edit costs between two graphs. The implementation is taken from
[FORI-GED](https://github.com/meffertj/FORI-GED).

For details on valid variable dimensions see [`validate_cost_function`](@ref), for
constructing the canonical cost function see [`get_default_edit_costs`](@ref).
"""
struct EditCosts
    c_ik::Array{Number,2}
    c_iε::Vector{Number}
    c_εk::Vector{Number}
    # note that for sparse graphs we might want to use a sparse array here
    c_ijkl::Array{Number,4}
    c_ijε::Array{Number,2}
    c_εkl::Array{Number,2}
end

"""
    validate_cost_function(c, G, H)

Validates that cost function `c` has appropriate dimensions for graphs `G` and `H`. This 
means

- `c.c_ik` has a cost for each mapping `i ∈ V(G)` to `k ∈ V(H)`
- `c.c_iε` has a cost for deleting each `i ∈ V(G)`
- `c.c_εk` has a cost for creating each `k ∈ V(H)`
- `c.c_ijkl` has a cost for each mapping from `ij ∈ E(G)` to  `kl ∈ E(H)`
- `c.c_ijε` has a cost for deleting each `ij ∈ E(G)`
- `c.c_εkl` has a cost for creating each `kl ∈ E(H)`
"""
function validate_cost_function(c::EditCosts, G::AbstractGraph, H::AbstractGraph)
    @assert size(c.c_ik) == (nv(G), nv(H))
    @assert size(c.c_iε) == (nv(G),)
    @assert size(c.c_εk) == (nv(H),)
    @assert size(c.c_ijkl) == (nv(G), nv(G), nv(H), nv(H))
    @assert size(c.c_ijε) == (nv(G), nv(G))
    @assert size(c.c_εkl) == (nv(H), nv(H))
    return nothing
end

"""
Returns the cost function where mappings are free, deleting or creating things costs
uniformly 1. This corresponds to "number of edits needed to go from `G` to `H`".
"""
function get_default_edit_costs(G::AbstractGraph, H::AbstractGraph)
    return EditCosts(
        zeros(Int, nv(G), nv(H)),
        ones(Int, nv(G)),
        ones(Int, nv(H)),
        zeros(Int, nv(G), nv(G), nv(H), nv(H)),
        ones(Int, nv(G), nv(G)),
        ones(Int, nv(H), nv(H)),
    )
end

"""
Adds objective function for F1 style formulations to `model`.
"""
function add_F1_objective!(
    model, c::EditCosts, vars::FullVariables, G::AbstractGraph, H::AbstractGraph
)
    @objective(
        model,
        Min,
        sum(vars.x[i, k] * c.c_ik[i, k] for i in 1:nv(G) for k in 1:nv(H)) +
            sum(vars.nodeDelG[i] * c.c_iε[i] for i in 1:nv(G)) +
            sum(vars.nodeDelH[k] * c.c_εk[k] for k in 1:nv(H)) +
            sum(
                vars.y[i, j, k, l] * c.c_ijkl[i, j, k, l] for i in 1:nv(G) for j in 1:nv(G)
                for k in 1:nv(H) for
                l in 1:nv(H) if has_edge(G, i, j) && has_edge(H, k, l) && i < j && k < l
            ) +
            sum(
                vars.edgeDelG[i, j] * c.c_ijε[i, j] for i in 1:nv(G) for
                j in 1:nv(G) if has_edge(G, i, j) && i < j
            ) +
            sum(
                vars.edgeDelH[k, l] * c.c_εkl[k, l] for k in 1:nv(H) for
                l in 1:nv(H) if has_edge(H, k, l) && k < l
            )
    )
    return nothing
end

"""
Adds objective function for F2 style formulations to `model`. `bidirectional` is used to
support input for [`FORI`](@ref) variables.
"""
function add_F2_objective!(
    model,
    c::EditCosts,
    vars::ReducedVariables,
    G::AbstractGraph,
    H::AbstractGraph,
    bidirectional::Bool=false,
)
    # Note that init=0 is necessary to allow for empty graph edge cases.
    K =
        sum(c.c_iε[i] for i in 1:nv(G); init=0) +
        sum(c.c_εk[k] for k in 1:nv(H); init=0) +
        sum(
            c.c_ijε[i, j] for i in 1:nv(G) for j in 1:nv(G) if has_edge(G, i, j) && i < j;
            init=0,
        ) +
        sum(
            c.c_εkl[k, l] for k in 1:nv(H) for l in 1:nv(H) if has_edge(H, k, l) && k < l;
            init=0,
        )
    @objective(
        model,
        Min,
        sum(
                vars.x[i, k] * (c.c_ik[i, k] - c.c_iε[i] - c.c_εk[k]) for i in 1:nv(G) for
                k in 1:nv(H);
                init=0,
            ) +
            sum(
                vars.y[i, j, k, l] * (c.c_ijkl[i, j, k, l] - c.c_ijε[i, j] - c.c_εkl[k, l])
                for i in 1:nv(G) for j in 1:nv(G) for k in 1:nv(H) for l in 1:nv(H) if
                has_edge(G, i, j) && has_edge(H, k, l) && i < j && (k < l || bidirectional);
                init=0,
            ) +
            K
    )
    return nothing
end

"""
Adds objective function for [`FORI`](@ref) formulation variables. The implementation is in
[`add_F2_objective!`](@ref).
"""
function add_FORI_objective!(
    model, c::EditCosts, vars::OrientedVariables, G::AbstractGraph, H::AbstractGraph
)
    return add_F2_objective!(model, c, ReducedVariables(vars.x, vars.z), G, H, true)
end

"""
Add node map constraints to model for F1 style formulations, i.e. each node is only mapped
to exactly one other node or deleted/created.
"""
function add_node_map_constraints!(model::Model, vars::FullVariables, G, H)
    @constraint(model, [i in 1:nv(G)], sum(vars.x[i, :]) + vars.nodeDelG[i] == 1)
    @constraint(model, [j in 1:nv(H)], sum(vars.x[:, j]) + vars.nodeDelH[j] == 1)
    return nothing
end

"""
Add node map constraints to model for F2 style and FORI formulations, i.e. each node is only
mapped to at most one other node (not being mapped implies being deleted or created).
"""
function add_node_map_constraints!(
    model::Model, vars::Union{ReducedVariables,OrientedVariables}, G, H
)
    @constraint(model, [i in 1:nv(G)], sum(vars.x[i, :]) <= 1)
    @constraint(model, [j in 1:nv(H)], sum(vars.x[:, j]) <= 1)
    return nothing
end

"""
Adds constraints on the edge map variables so each edge is mapped to exactly one other or be
deleted/created. Note that these constraints are only necessary in F1 style formulations, as
they are implied in F2 styled and FORI.
"""
function add_edge_map_constraints!(model::Model, vars::FullVariables, G, H)
    # Note that init=0 is necessary to allow for empty graph edge cases.
    @constraint(
        model,
        [i in vertices(G), j in vertices(G); has_edge(G, i, j) && i < j],
        sum(vars.y[i, j, :, :]; init=0) + vars.edgeDelG[i, j] == 1
    )
    @constraint(
        model,
        [k in vertices(H), l in vertices(H); has_edge(H, k, l) && k < l],
        sum(vars.y[:, :, k, l]; init=0) + vars.edgeDelH[k, l] == 1
    )
    return nothing
end

"""
Add simplest form of topology constraints: If edge `ij` is mapped to `kl`, then `i` or `j`
must map to `k`, and `i` or `j` must map to `l`.
"""
function add_simple_topology_constraints!(model::Model, vars::Variables, G, H)
    @constraint(
        model,
        [
            i in vertices(G),
            j in vertices(G),
            k in vertices(H),
            l in vertices(H);
            has_edge(G, i, j) && has_edge(H, k, l) && i < j && k < l,
        ],
        vars.y[i, j, k, l] <= vars.x[i, k] + vars.x[j, k]
    )
    @constraint(
        model,
        [
            i in vertices(G),
            j in vertices(G),
            k in vertices(H),
            l in vertices(H);
            has_edge(G, i, j) && has_edge(H, k, l) && i < j && k < l,
        ],
        vars.y[i, j, k, l] <= vars.x[i, l] + vars.x[j, l]
    )
    return nothing
end

"""
Add improved topology constraints, replacing [`add_simple_topology_constraints!`](@ref): 
If edge `ij` is mapped to any edge incident to `k`, then `i` or `j` must be mapped to `k`.
"""
function add_improved_topology_constraints_G_to_H!(model::Model, vars::Variables, G, H)
    @constraint(
        model,
        [i in vertices(G), j in vertices(G), k in vertices(H); has_edge(G, i, j) && i < j],
        sum(vars.y[i, j, k, :]; init=0) + sum(vars.y[i, j, :, k]; init=0) <=
            vars.x[i, k] + vars.x[j, k]
    )
    return nothing
end

"""
Add topology constraints mirroring [`add_improved_topology_constraints_G_to_H!`](@ref), but
backwards: If any edge incident to `i` is mapped to `kl`, then `i` must be mapped to `k` or
`l`.
"""
function add_improved_topology_constraints_H_to_G!(model::Model, vars::Variables, G, H)
    @constraint(
        model,
        [k in vertices(H), l in vertices(H), i in vertices(G); has_edge(H, k, l) && k < l],
        sum(vars.y[i, :, k, l]; init=0) + sum(vars.y[:, i, k, l]; init=0) <=
            vars.x[i, k] + vars.x[i, l]
    )
    return nothing
end

"""
Add topology constraints for FORI. Uses the same implications as used in
[`add_improved_topology_constraints_G_to_H!`](@ref) and
[`add_improved_topology_constraints_H_to_G!`](@ref), but since edges in `G` are oriented and
each edge in `H` has two possible orientations, the implications are more precise.
"""
function add_oriented_topology_constraints!(model::Model, vars::OrientedVariables, G, H)
    @constraint(
        model,
        [k in vertices(H), i in vertices(G), j in vertices(G); has_edge(G, i, j) && i < j],
        sum(vars.z[i, j, k, :]; init=0) <= vars.x[i, k]
    )
    @constraint(
        model,
        [k in vertices(H), i in vertices(G), j in vertices(G); has_edge(G, i, j) && i < j],
        sum(vars.z[i, j, :, k]; init=0) <= vars.x[j, k]
    )
    @constraint(
        model,
        [k in vertices(H), l in vertices(H), i in vertices(G); has_edge(H, k, l)],
        sum(vars.z[i, :, k, l]; init=0) + sum(vars.z[:, i, l, k]; init=0) <= vars.x[i, k]
    )
    return nothing
end

"""
    construct_formulation!(::Type{<:Formulation}, model, G, H;
        c)

Modify `model` to use a specific formulation to solve the graph edit distance problem. Each
formulation has its own method implementing the required variables, constraints and
objective.
"""
function construct_formulation! end

function construct_formulation!(
    ::Type{F1}, model, G, H, c::EditCosts=get_default_edit_costs(G, H)
)
    vars = create_model_vars_full!(model, G, H)

    add_node_map_constraints!(model, vars, G, H)
    add_edge_map_constraints!(model, vars, G, H)

    add_simple_topology_constraints!(model, vars, G, H)

    add_F1_objective!(model, c, vars, G, H)
    return nothing
end

function construct_formulation!(
    ::Type{F2minus}, model, G, H, c::EditCosts=get_default_edit_costs(G, H)
)
    vars = create_model_vars_reduced!(model, G, H)

    add_node_map_constraints!(model, vars, G, H)
    # Edgemap constraints are implied, so we can skip them.

    add_simple_topology_constraints!(model, vars, G, H)

    add_F2_objective!(model, c, vars, G, H)
    return nothing
end

function construct_formulation!(
    ::Type{F2}, model, G, H, c::EditCosts=get_default_edit_costs(G, H)
)
    vars = create_model_vars_reduced!(model, G, H)

    add_node_map_constraints!(model, vars, G, H)
    # Edgemap constraints are implied, so we can skip them.

    add_improved_topology_constraints_G_to_H!(model, vars, G, H)

    add_F2_objective!(model, c, vars, G, H)
    return nothing
end

function construct_formulation!(
    ::Type{F2plus}, model, G, H, c::EditCosts=get_default_edit_costs(G, H)
)
    vars = create_model_vars_reduced!(model, G, H)

    add_node_map_constraints!(model, vars, G, H)
    # Edgemap constraints are implied, so we can skip them.

    add_improved_topology_constraints_G_to_H!(model, vars, G, H)
    add_improved_topology_constraints_H_to_G!(model, vars, G, H)

    add_F2_objective!(model, c, vars, G, H)
    return nothing
end

function construct_formulation!(
    ::Type{F1prime}, model, G, H, c::EditCosts=get_default_edit_costs(G, H)
)
    vars = create_model_vars_full!(model, G, H)

    add_node_map_constraints!(model, vars, G, H)
    add_edge_map_constraints!(model, vars, G, H)

    add_improved_topology_constraints_G_to_H!(model, vars, G, H)

    add_F1_objective!(model, c, vars, G, H)
    return nothing
end

function construct_formulation!(
    ::Type{F1plus}, model, G, H, c::EditCosts=get_default_edit_costs(G, H)
)
    vars = create_model_vars_full!(model, G, H)

    add_node_map_constraints!(model, vars, G, H)
    add_edge_map_constraints!(model, vars, G, H)

    add_improved_topology_constraints_G_to_H!(model, vars, G, H)
    add_improved_topology_constraints_H_to_G!(model, vars, G, H)

    add_F1_objective!(model, c, vars, G, H)
    return nothing
end

function construct_formulation!(
    ::Type{FORI}, model, G, H, c::EditCosts=get_default_edit_costs(G, H)
)
    vars = create_model_vars_bidirectional!(model, G, H)
    add_node_map_constraints!(model, vars, G, H)
    # Edgemap constraints are implied, so we can skip them.

    add_oriented_topology_constraints!(model, vars, G, H)

    add_FORI_objective!(model, c, vars, G, H)
    return nothing
end

"""
    edit_distance!(model, G, H;
        c, formulation
    )

Modify a JuMP model to compute the graph edit distance between undirected graphs `G` and `H`
given cost function `c` using formulation `formulation`. See [`edit_distance`](@ref) for
more details.
"""
function edit_distance!(
    model,
    G::SimpleGraph,
    H::SimpleGraph;
    c::EditCosts=get_default_edit_costs(G, H),
    formulation::Type{<:Formulation}=FORI,
)
    construct_formulation!(formulation, model, G, H, c)
    if num_variables(model) == 0
        # in this case, one or both of the graphs have no edges and/or no nodes.
        # For the solver to work, we add an unbounded dummy variable
        @variable(model, dummy)
    end
    return nothing
end

"""
    edit_distance(G, H;
        c, formulation, optimizer
    )

Compute the graph edit distance between undirected graphs `G` and `H` given edit costs `c`.

# Arguments
- `G::Graphs.SimpleGraph`
- `H::Graphs.SimpleGraph`

# Keywords
- `c::EditCosts`: pairwise edit cost struct for the graphs. See
    [`GraphsOptim.EditCosts`](@ref). Defaults to the canonical edit costs
    [`GraphsOptim.get_default_edit_costs`](@ref).
- `formulation::Type{<:Formulation}`: the ILP formulation to use. Defaults to the most 
    powerful `FORI`.
- `optimizer`: JuMP-compatible solver (default is `HiGHS.Optimizer`)

# Returns
A named tuple containing:

- `objective_value`: the optimal graph edit distance.
- `node_matching`: the vertex-matching matrix associated with the optimal solution.
"""
function edit_distance(
    G::SimpleGraph,
    H::SimpleGraph;
    c::EditCosts=get_default_edit_costs(G, H),
    formulation::Type{<:Formulation}=FORI,
    optimizer=HiGHS.Optimizer,
)
    model = Model(optimizer)
    set_silent(model)
    edit_distance!(model, G, H; c=c, formulation=formulation)
    optimize!(model)
    if termination_status(model) != OPTIMAL
        error("Graph edit distance was not solved optimally.")
    end
    return (; objective_value=objective_value(model), node_matching=value.(model[:x]))
end
