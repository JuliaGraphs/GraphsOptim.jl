using GraphsOptim
using Graphs
using Test

g = Graphs.random_regular_graph(10, 5)

for _ in 1:10
    vertex_weights = rand(nv(g))
    stable = GraphsOptim.maximum_weight_independent_set(g; vertex_weights=vertex_weights)
    if length(stable) > 1
        for idx in 1:(length(stable) - 1)
            @test !Graphs.has_edge(g, stable[idx], stable[idx + 1])
        end
    end
end

g2 = complete_graph(3)
add_vertex!(g2)
add_edge!(g2, 3, 4)
stable = GraphsOptim.maximum_weight_independent_set(g2)
@test length(stable) == 2
@test 4 in stable
