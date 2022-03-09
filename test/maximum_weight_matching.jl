@testset "Maximum weight matching" begin
    g = complete_graph(3)
    w = [
        1 2 1
        1 1 1
        3 1 1
    ]
    match = maximum_weight_matching(g, w; optimizer=GLPK.Optimizer)
    @test match.mate[1] == 3
    @test match.weight ≈ 3

    g = complete_graph(3)
    w = zeros(3, 3)
    w[1, 2] = 1
    w[3, 2] = 1
    w[1, 3] = 1
    match = maximum_weight_matching(g, w; optimizer=GLPK.Optimizer)
    @test match.weight ≈ 1

    g = Graph(4)
    add_edge!(g, 1, 3)
    add_edge!(g, 1, 4)
    add_edge!(g, 2, 4)

    w = zeros(4, 4)
    w[1, 3] = 1
    w[1, 4] = 3
    w[2, 4] = 1

    match = maximum_weight_matching(g, w; optimizer=GLPK.Optimizer)
    @test match.weight ≈ 3
    @test match.mate[1] == 4
    @test match.mate[2] == -1
    @test match.mate[3] == -1
    @test match.mate[4] == 1

    g = Graph(4)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 1)
    add_edge!(g, 3, 4)
    match = maximum_weight_matching(g; optimizer=GLPK.Optimizer)
    @test match.weight ≈ 2
    @test match.mate[1] == 2
    @test match.mate[2] == 1
    @test match.mate[3] == 4
    @test match.mate[4] == 3

    w = zeros(4, 4)
    w[1, 2] = 1
    w[2, 3] = 1
    w[1, 3] = 1
    w[3, 4] = 1

    match = maximum_weight_matching(g, w; optimizer=GLPK.Optimizer)
    @test match.weight ≈ 2
    @test match.mate[1] == 2
    @test match.mate[2] == 1
    @test match.mate[3] == 4
    @test match.mate[4] == 3

    w = zeros(4, 4)
    w[1, 2] = 1
    w[2, 3] = 1
    w[1, 3] = 5
    w[3, 4] = 1

    match = maximum_weight_matching(g, w; optimizer=GLPK.Optimizer)
    @test match.weight ≈ 5
    @test match.mate[1] == 3
    @test match.mate[2] == -1
    @test match.mate[3] == 1
    @test match.mate[4] == -1
end

@testset "Maximum weight maximal matching" begin
    g = complete_bipartite_graph(2, 2)
    w = zeros(4, 4)
    w[1, 3] = 10.0
    w[1, 4] = 1.0
    w[2, 3] = 2.0
    w[2, 4] = 11.0
    match = maximum_weight_maximal_matching(g, w; optimizer=GLPK.Optimizer)
    @test match.weight ≈ 21
    @test match.mate[1] == 3
    @test match.mate[3] == 1
    @test match.mate[2] == 4
    @test match.mate[4] == 2

    g = complete_bipartite_graph(2, 4)
    w = zeros(6, 6)
    w[1, 3] = 10
    w[1, 4] = 0.5
    w[2, 3] = 11
    w[2, 4] = 1
    match = maximum_weight_maximal_matching(g, w; optimizer=GLPK.Optimizer)
    @test match.weight ≈ 11.5
    @test match.mate[1] == 4
    @test match.mate[4] == 1
    @test match.mate[2] == 3
    @test match.mate[3] == 2

    g = complete_bipartite_graph(2, 6)
    w = zeros(8, 8)
    w[1, 3] = 10
    w[1, 4] = 0.5
    w[2, 3] = 11
    w[2, 4] = 1
    w[2, 5] = -1
    w[2, 6] = -1
    match = maximum_weight_maximal_matching(g, w, 0; optimizer=GLPK.Optimizer)
    @test match.weight ≈ 11.5
    @test match.mate[1] == 4
    @test match.mate[4] == 1
    @test match.mate[2] == 3
    @test match.mate[3] == 2

    g = complete_bipartite_graph(4, 2)
    w = zeros(6, 6)
    w[3, 5] = 10
    w[3, 6] = 0.5
    w[2, 5] = 11
    w[1, 6] = 1
    w[1, 5] = -1

    match = maximum_weight_maximal_matching(g, w, 0; optimizer=GLPK.Optimizer)
    @test match.weight ≈ 12
    @test match.mate[1] == 6
    @test match.mate[2] == 5
    @test match.mate[3] == -1
    @test match.mate[4] == -1
    @test match.mate[5] == 2
    @test match.mate[6] == 1
end
