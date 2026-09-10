using GraphsOptim, Graphs, JuMP, HiGHS
using Test

all_formulations = [F1, F1prime, F1plus, F2, F2minus, F2plus, FORI]

@testset "Small Graph" begin
    G = Graph(3)
    add_edge!(G, 1, 2)

    H = Graph(3)
    add_edge!(H, 2, 3)
    add_edge!(H, 1, 3)

    for formulation in all_formulations
        model = Model(HiGHS.Optimizer)
        GraphsOptim.edit_distance!(model, G, H; formulation=formulation)
        set_silent(model)
        optimize!(model)
        @test objective_value(model) == 1
        # testing the normal entrypoint
        @test GraphsOptim.edit_distance(G, H).objective_value == 1
    end
end

@testset "Cost function validations" begin
    G = Graph(3)
    add_edge!(G, 1, 2)

    H = Graph(3)
    add_edge!(H, 2, 3)
    add_edge!(H, 1, 3)

    c = GraphsOptim.get_default_edit_costs(G, H)
    @test isnothing(GraphsOptim.validate_cost_function(c, G, H))

    broken_c = GraphsOptim.EditCosts(
        zeros(Int, 1, 2), c.c_iε, c.c_εk, c.c_ijkl, c.c_ijε, c.c_εkl
    )
    @test_throws AssertionError GraphsOptim.validate_cost_function(broken_c, G, H)

    Hprime = Graph(4)
    @test_throws AssertionError GraphsOptim.validate_cost_function(c, G, Hprime)
end

@testset "Cost Functions Simple" begin
    G = Graph(3)
    add_edge!(G, 1, 2)

    H = Graph(3)
    add_edge!(H, 2, 3)
    add_edge!(H, 1, 3)
    c_origin = GraphsOptim.get_default_edit_costs(G, H)

    # increating costs for node subsitutions 
    c = GraphsOptim.EditCosts(
        ones(Int, nv(G), nv(H)),
        5 * c_origin.c_iε,
        5 * c_origin.c_εk,
        c_origin.c_ijkl,
        c_origin.c_ijε,
        c_origin.c_εkl,
    )
    for formulation in all_formulations
        model = Model(HiGHS.Optimizer)
        GraphsOptim.edit_distance!(model, G, H; c=c, formulation=formulation)
        set_silent(model)
        optimize!(model)
        @test objective_value(model) == 4
    end
    # increating costs for edge subsitutions 
    c = GraphsOptim.EditCosts(
        c_origin.c_ik,
        c_origin.c_iε,
        c_origin.c_εk,
        c_origin.c_ijkl,
        100 * c_origin.c_ijε,
        10 * c_origin.c_εkl,
    )
    for formulation in all_formulations
        model = Model(HiGHS.Optimizer)
        GraphsOptim.edit_distance!(model, G, H; c=c, formulation=formulation)
        set_silent(model)
        optimize!(model)
        @test objective_value(model) == 10
    end

    # force a specific suboptimal node map with edit costs 
    G = Graph(3)
    add_edge!(G, 2, 3)

    H = Graph(4)
    add_edge!(H, 1, 3)
    add_edge!(H, 3, 4)
    c_origin = GraphsOptim.get_default_edit_costs(G, H)
    custom_substitution_cost = 100 .+ c_origin.c_ik
    custom_substitution_cost[1, 1] = 0
    custom_substitution_cost[2, 2] = 0
    custom_substitution_cost[3, 3] = 0
    c = GraphsOptim.EditCosts(
        custom_substitution_cost,
        c_origin.c_iε,
        c_origin.c_εk,
        c_origin.c_ijkl,
        c_origin.c_ijε,
        c_origin.c_εkl,
    )
    for formulation in all_formulations
        model = Model(HiGHS.Optimizer)
        vars = GraphsOptim.edit_distance!(model, G, H; c=c, formulation=formulation)
        set_silent(model)
        optimize!(model)
        @test objective_value(model) == 4
        @test value(model[:x][1, 1]) == 1
        @test value(model[:x][2, 2]) == 1
        @test value(model[:x][3, 3]) == 1
    end
end

@testset "NodeMatchingOutput" begin
    G = path_graph(3)
    H = path_graph(4)
    dist, node_matching = GraphsOptim.edit_distance(G, H)
    @test dist == 2

    # test the matrix dimensions and that the rounded result is a 0-1 matrix
    @test size(node_matching) == (3, 4)
    rounded_node_matching = round.(node_matching)
    @test all(x -> x in (0, 1), rounded_node_matching)

    # test its a matching matrix, accounting for the fact that it isn't a perfect match
    for i in 1:3
        @test sum(rounded_node_matching[i, :]) == 1
    end
    @test sum(rounded_node_matching[:, 1] + rounded_node_matching[:, 4]) == 1
    for j in 2:3
        @test sum(rounded_node_matching[:, j]) == 1
    end
end

@testset "Edgecases" begin
    G = Graph(3)

    H = Graph(3)
    add_edge!(H, 1, 2)

    # one graph has no edges
    for formulation in all_formulations
        model = Model(HiGHS.Optimizer)
        GraphsOptim.edit_distance!(model, G, H; formulation=formulation)
        set_silent(model)
        optimize!(model)
        @test objective_value(model) == 1

        # reverse order
        model = Model(HiGHS.Optimizer)
        GraphsOptim.edit_distance!(model, H, G; formulation=formulation)
        set_silent(model)
        optimize!(model)
        @test objective_value(model) == 1
    end

    # both graphs have no edges
    H = Graph(5)
    for formulation in all_formulations
        model = Model(HiGHS.Optimizer)
        GraphsOptim.edit_distance!(model, G, H; formulation=formulation)
        set_silent(model)
        optimize!(model)
        @test objective_value(model) == 2
    end

    # no nodes in one graph
    for formulation in all_formulations
        model = Model(HiGHS.Optimizer)
        GraphsOptim.edit_distance!(model, G, Graph(0); formulation=formulation)
        set_silent(model)
        optimize!(model)
        @test objective_value(model) == 3
    end
    # no nodes in both graphs
    for formulation in all_formulations
        model = Model(HiGHS.Optimizer)
        GraphsOptim.edit_distance!(model, Graph(0), Graph(0); formulation=formulation)
        set_silent(model)
        optimize!(model)
        @test objective_value(model) == 0

        # the same behaviour in the edit_distance function
        @test GraphsOptim.edit_distance(Graph(0), Graph(0)).objective_value == 0
    end
end

@testset "Invalid Input" begin
    G = Graph(3)
    H = DiGraph(4)
    @test_throws MethodError GraphsOptim.edit_distance(G, H)
end
