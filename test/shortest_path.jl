using Graphs
using GraphsOptim
using Test

num_vertices = 20

function test_shortest_path(g)
    edge_cost = rand(1:100, num_vertices, num_vertices) .* adjacency_matrix(g)

    source = rand(1:num_vertices)
    target = rand(setdiff(Set(1:num_vertices), Set(source)))

    dijkstra_parents = Graphs.dijkstra_shortest_paths(g, source, edge_cost).parents
    dijkstra_ans = [target]
    while true
        parent_vx = dijkstra_parents[dijkstra_ans[begin]]
        pushfirst!(dijkstra_ans, parent_vx)
        if parent_vx == source
            break
        end
    end

    mp_ans = GraphsOptim.shortest_path(g, source, target, edge_cost)

    @test all(dijkstra_ans .== mp_ans)
end

g = Graphs.SimpleDiGraph(num_vertices, 5 * num_vertices; seed=0)
test_shortest_path(g)

g = Graphs.SimpleGraph(num_vertices, 3 * num_vertices; seed=0)
test_shortest_path(g)
