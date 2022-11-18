g = complete_graph(3)
w = [
    1 2 1
    1 1 1
    3 1 1
]
match = maximum_weight_matching(g, w; optimizer=HiGHS.Optimizer)
@test match.mate[1] == 3
@test match.weight ≈ 3

g = complete_graph(3)
w = zeros(3, 3)
w[1, 2] = 1
w[3, 2] = 1
w[1, 3] = 1
match = maximum_weight_matching(g, w; optimizer=HiGHS.Optimizer)
@test match.weight ≈ 1

g = Graph(4)
add_edge!(g, 1, 3)
add_edge!(g, 1, 4)
add_edge!(g, 2, 4)

w = zeros(4, 4)
w[1, 3] = 1
w[1, 4] = 3
w[2, 4] = 1

match = maximum_weight_matching(g, w; optimizer=HiGHS.Optimizer)
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
match = maximum_weight_matching(g; optimizer=HiGHS.Optimizer)
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

match = maximum_weight_matching(g, w; optimizer=HiGHS.Optimizer)
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

match = maximum_weight_matching(g, w; optimizer=HiGHS.Optimizer)
@test match.weight ≈ 5
@test match.mate[1] == 3
@test match.mate[2] == -1
@test match.mate[3] == 1
@test match.mate[4] == -1
