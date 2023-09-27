using Graphs
using GraphsOptim
using Test

num_vertices = 20
g = Graphs.SimpleDiGraph(num_vertices, 5 * num_vertices; seed = 0)
edge_cost = rand(1:100, num_vertices, num_vertices) .* adjacency_matrix(g)

src = rand(1:num_vertices)
dst = rand(setdiff(Set(1:num_vertices), Set(src)))

dijkstra_parents = Graphs.dijkstra_shortest_paths(g, src, edge_cost).parents
dijkstra_ans = [dst]
while true
	parent_vx = dijkstra_parents[dijkstra_ans[begin]]
	pushfirst!(dijkstra_ans, parent_vx)
	if parent_vx == src
		break
	end
end

mp_ans = GraphsOptim.shortest_path(g, src, dst, edge_cost)

@test all(dijkstra_ans .== mp_ans)