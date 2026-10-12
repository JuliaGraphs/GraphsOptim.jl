using Graphs
using GraphsOptim
using SimpleWeightedGraphs
using Test

function is_isomorphism(g, h, f)
    return Set(keys(f)) == Set(vertices(g)) &&
           Set(values(f)) == Set(vertices(h)) &&
           all(
               has_edge(g, u, v) == has_edge(h, f[u], f[v]) for u in vertices(g) for
               v in vertices(g)
           )
end

@testset "Undirected graphs" begin
    g = path_graph(5)
    h = SimpleGraph(5)
    permutation = [3, 5, 1, 4, 2]
    for edge in edges(g)
        add_edge!(h, permutation[src(edge)], permutation[dst(edge)])
    end

    f = graph_isomorphim(g, h)
    @test is_isomorphism(g, h, f)
    @test graph_isomorphim(SimpleGraph(0), SimpleGraph(0)) == Dict{Int,Int}()

    two_triangles = SimpleGraph(6)
    for triangle in ([1, 2, 3], [4, 5, 6])
        add_edge!(two_triangles, triangle[1], triangle[2])
        add_edge!(two_triangles, triangle[2], triangle[3])
        add_edge!(two_triangles, triangle[3], triangle[1])
    end
    @test_throws ArgumentError graph_isomorphim(cycle_graph(6), two_triangles)
    @test_throws ArgumentError graph_isomorphim(path_graph(3), complete_graph(3))
end

@testset "Directed graphs" begin
    g = SimpleDiGraph(4)
    add_edge!.(Ref(g), [1, 1, 2, 3], [2, 3, 4, 4])
    h = SimpleDiGraph(4)
    permutation = [4, 2, 1, 3]
    for edge in edges(g)
        add_edge!(h, permutation[src(edge)], permutation[dst(edge)])
    end

    @test is_isomorphism(g, h, graph_isomorphim(g, h))
    @test_throws ArgumentError graph_isomorphim(g, SimpleGraph(g))
end

@testset "Weighted graphs" begin
    A = [0.0 1.5 0.0 0.0; 1.5 0.0 2.5 0.0; 0.0 2.5 0.0 4.5; 0.0 0.0 4.5 0.0]
    g = SimpleWeightedGraph(A)
    permutation = [3, 1, 4, 2]
    B = zeros(4, 4)
    B[permutation, permutation] = A
    h = SimpleWeightedGraph(B)

    f = graph_isomorphim(g, h)
    @test all(
        weights(g)[u, v] == weights(h)[f[u], f[v]] for u in vertices(g) for v in vertices(g)
    )

    B[permutation[1], permutation[2]] = 10.0
    B[permutation[2], permutation[1]] = 10.0
    @test_throws ArgumentError graph_isomorphim(g, SimpleWeightedGraph(B))
end
