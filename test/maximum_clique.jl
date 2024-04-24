using GraphsOptim
using Graphs
using Test

g = Graphs.random_regular_graph(10, 5)

for _ in 1:10
    vertex_weights = rand(nv(g))
    clique = GraphsOptim.maximum_weight_clique(g; vertex_weights=vertex_weights)
    if length(clique) > 1
        for idx in 1:(length(clique) - 1)
            @test Graphs.has_edge(g, clique[idx], clique[idx + 1])
        end
    end
end

g2 = complete_graph(3)
add_vertex!(g2)
add_edge!(g2, 3, 4)
clique = GraphsOptim.maximum_weight_clique(g2)
@test sort(clique) == 1:3
