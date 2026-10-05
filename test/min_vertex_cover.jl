using Graphs: SimpleGraph
using GraphsOptim
using Test

adj_matrix = [
    0 1 1 0 0 0
    1 0 1 1 1 1
    1 1 0 0 0 0
    0 1 0 0 0 0
    0 1 0 0 0 0
    0 1 0 0 0 0
]

expected_ans_options = [Set([1, 2]), Set([2, 3])]

graph = SimpleGraph(adj_matrix)

@test Set(min_vertex_cover(graph)) in expected_ans_options

# Covers must be valid even when the solver returns values slightly below 1
using Graphs: erdos_renyi, edges, src, dst, nv
for seed in 1:5
    g = erdos_renyi(50, 0.1; seed)
    cover = Set(min_vertex_cover(g))
    @test all(src(e) in cover || dst(e) in cover for e in edges(g))
    @test length(cover) == nv(g) - length(maximum_weight_independent_set(g))
end
