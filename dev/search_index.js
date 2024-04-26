var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"EditURL = \"https://github.com/JuliaGraphs/GraphsOptim.jl/blob/main/README.md\"","category":"page"},{"location":"#GraphsOptim","page":"Home","title":"GraphsOptim","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Dev) (Image: Build Status) (Image: Coverage) (Image: Code Style: Blue)","category":"page"},{"location":"","page":"Home","title":"Home","text":"A package for graph optimization algorithms that rely on mathematical programming.","category":"page"},{"location":"#Getting-started","page":"Home","title":"Getting started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package is still experimental, which is why it is not yet registered. To install it, you need to use the GitHub URL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg; Pkg.add(url=\"https://github.com/JuliaGraphs/GraphsOptim.jl\")","category":"page"},{"location":"#Roadmap","page":"Home","title":"Roadmap","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package only contains a few algorithms, and we would like to add more. New contributors are always welcome: just pick a problem from our roadmap issue, open a pull request following the guidelines, and we will help you get it merged!","category":"page"},{"location":"algorithms/#Algorithms","page":"Algorithms","title":"Algorithms","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"GraphsOptim","category":"page"},{"location":"algorithms/#GraphsOptim","page":"Algorithms","title":"GraphsOptim","text":"GraphsOptim\n\nA package for graph optimization algorithms that rely on mathematical programming.\n\n\n\n\n\n","category":"module"},{"location":"algorithms/#Flow","page":"Algorithms","title":"Flow","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"min_cost_flow\nGraphsOptim.min_cost_flow!","category":"page"},{"location":"algorithms/#GraphsOptim.min_cost_flow","page":"Algorithms","title":"GraphsOptim.min_cost_flow","text":"min_cost_flow(\n    g, vertex_demand, edge_cost, edge_min_capacity, edge_max_capacity;\n    integer, optimizer\n)\n\nCompute a minimum cost flow over a directed graph.\n\nArguments\n\ng::Graphs.AbstractGraph: a directed graph G = (V, E)\nvertex_demand::AbstractVector: a vector in Rⱽ giving the flow requested by each vertex (should be positive for sinks, negative for sources and zero elsewhere)\nedge_cost::AbstractMatrix: a vector in Rᴱ giving the cost of a unit of flow on each edge\nedge_min_capacity::AbstractMatrix: a vector in Rᴱ giving the minimum flow allowed on each edge\nedge_max_capacity::AbstractMatrix: a vector in Rᴱ giving the maximum flow allowed on each edge\n\nKeyword arguments\n\ninteger::Bool: whether the flow should be integer-valued or real-valued\noptimizer: JuMP-compatible solver (default is HiGHS.Optimizer)\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.min_cost_flow!","page":"Algorithms","title":"GraphsOptim.min_cost_flow!","text":"min_cost_flow!(\n    model,\n    g, vertex_demand, edge_cost, edge_min_capacity, edge_max_capacity;\n    var_name, integer\n)\n\nModify a JuMP model by adding the variable, constraints and objective necessary to compute a minimum cost flow over a directed graph.\n\nThe flow variable will be named var_name, see min_cost_flow for details on the other arguments.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"We denote by:","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"f the edge flow variable\nc the edge cost\na and b the min and max edge capacity\nd the vertex demand","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"The objective function is","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"min_f in mathbbR^E sum_(u v) in E c(u v) f(u v)","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"The edge capacity constraint dictates that for all (u v) in E,","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"a(u v) leq f(u v) leq b(u v)","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"The flow conservation constraint with node demand dictates that for all v in V,","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"f^-(v) = d(v) + f^+(v)","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"where the incoming flow f^-(v) and outgoing flow f^+(v) are defined as","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"f^-(v) = sum_u in N^-(v) f(u v) quad textand quad f^+(v) = sum_w in N^+(v) f(v w)","category":"page"},{"location":"algorithms/#Shortest-Path","page":"Algorithms","title":"Shortest Path","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"shortest_path\nGraphsOptim.shortest_path!","category":"page"},{"location":"algorithms/#GraphsOptim.shortest_path","page":"Algorithms","title":"GraphsOptim.shortest_path","text":"shortest_path(\n\tg, source, target, edge_cost;\n    integer, optimizer\n)\n\nCompute the shortest path from a source to a target vertex over a graph. Returns a sequence of vertices from source to destination.\n\nArguments\n\ng::Graphs.AbstractGraph: a directed or undirected graph G = (V, E)\nsource::Int: source vertex of the path\ntarget::Int: target vertex of the path\nedge_cost::AbstractMatrix: a vector in Rᴱ giving the cost of traversing the edge\n\nKeyword arguments\n\ninteger::Bool: whether the path should be integer-valued or real-valued (default is true)\noptimizer: JuMP-compatible solver (default is HiGHS.Optimizer)\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.shortest_path!","page":"Algorithms","title":"GraphsOptim.shortest_path!","text":"shortest_path!(\n\tmodel,\n\tg, source, target, edge_cost;\n    var_name, integer\n)\n\nModify a JuMP model by adding the variable, constraints and objective necessary to compute the shortest path between 2 vertices over a graph.\n\nThe edge selection variable will be named var_name, see shortest_path for details on the other arguments.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"A special case of minimum cost flow without edge capacities, and where vertex demands are 0 everywhere except at the source (-1) and target (+1).","category":"page"},{"location":"algorithms/#Assignment","page":"Algorithms","title":"Assignment","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"danger: Work in progress\nCome back later!","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"min_cost_assignment\nGraphsOptim.min_cost_assignment!","category":"page"},{"location":"algorithms/#GraphsOptim.min_cost_assignment","page":"Algorithms","title":"GraphsOptim.min_cost_assignment","text":"min_cost_assignment(\n    edge_cost;\n    optimizer\n)\n\nCompute a minimum cost assignment over a bipartite graph.\n\nArguments\n\nedge_cost::AbstractMatrix: a matrix in Rᵁˣⱽ giving the cost of matching u ∈ U to v ∈ V\n\nKeyword arguments\n\ninteger::Bool: whether the flow should be integer-valued or real-valued\noptimizer: JuMP-compatible solver (default is HiGHS.Optimizer)\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.min_cost_assignment!","page":"Algorithms","title":"GraphsOptim.min_cost_assignment!","text":"min_cost_assignment!(\n    model,\n    edge_cost;\n    var_name, integer\n)\n\nModify a JuMP model by adding the variable, constraints and objective necessary to compute a minimum cost assignment over a bipartite graph.\n\nThe assignment variable will be named var_name, see min_cost_assignment for details on the other arguments.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#Minimum-Vertex-Cover","page":"Algorithms","title":"Minimum Vertex Cover","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"min_vertex_cover\nGraphsOptim.min_vertex_cover!","category":"page"},{"location":"algorithms/#GraphsOptim.min_vertex_cover","page":"Algorithms","title":"GraphsOptim.min_vertex_cover","text":"min_vertex_cover(\n    g, optimizer\n)\n\nCompute a minimum vertex cover of an undirected graph.\n\nArguments\n\ng::Graphs.AbstractGraph: an undirected graph G = (V, E)\n\nKeyword arguments\n\noptimizer: JuMP-compatible solver (default is HiGHS.Optimizer)\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.min_vertex_cover!","page":"Algorithms","title":"GraphsOptim.min_vertex_cover!","text":"min_vertex_cover!(\n    model, g;\n    var_name\n)\n\nModify a JuMP model by adding the variable, constraints and objective necessary to compute a minimum vertex cove of an undirected graph.\n\nThe vertex cover indicator variable will be named var_name\n\n\n\n\n\n","category":"function"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"Finds a subset S subset V of vertices of an undirected graph G = (VE) such that forall (uv) in E u in S lor v in S","category":"page"},{"location":"algorithms/#Maximum-Weight-Independent-Set","page":"Algorithms","title":"Maximum Weight Independent Set","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"maximum_weight_independent_set\nGraphsOptim.maximum_weight_independent_set!","category":"page"},{"location":"algorithms/#GraphsOptim.maximum_weight_independent_set","page":"Algorithms","title":"GraphsOptim.maximum_weight_independent_set","text":"maximum_weight_independent_set(g; optimizer, binary, vertex_weights)\n\nComputes a maximum-weighted independent set or stable set of g.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.maximum_weight_independent_set!","page":"Algorithms","title":"GraphsOptim.maximum_weight_independent_set!","text":"maximum_independent_set!(model, g; var_name)\n\nComputes in-place in the JuMP model a maximum-weighted independent set of g. An optional vertex_weights vector can be passed to the graph, defaulting to uniform weights (computing a maximum size independent set).\n\n\n\n\n\n","category":"function"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"Finds a subset S subset V of vertices of maximal weight of an undirected graph G = (VE) such that forall (uv) in E u notin S lor v notin S.","category":"page"},{"location":"algorithms/#Graph-matching","page":"Algorithms","title":"Graph matching","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"danger: Work in progress\nCome back later!","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"graph_matching\nGraphsOptim.graph_matching_step_size","category":"page"},{"location":"algorithms/#GraphsOptim.graph_matching","page":"Algorithms","title":"GraphsOptim.graph_matching","text":"graph_matching(\n    algo, A, B;\n    optimizer, P_init, max_iter, tol,\n    max_iter_sinkhorn, regularizer\n)\n\nCompute an approximately optimal alignment between two graphs using one of two variants of Frank-Wolfe (FAQ or GOAT)\n\nThe output is a tuple that contains:\n\nthe permutation matrix P defining the alignment;\nthe distance between the permuted graphs;\na boolean indicating if the algorithm converged.\n\nArguments\n\nalgo: allows dispatch based on the choice of algorithm, either FAQ() or GOAT()\nA: the adjacency matrix of the first graph\nB: the adjacency matrix of the second graph\n\nKeyword arguments\n\nFor both algorithms:\n\noptimizer: JuMP-compatible solver (default is HiGHS.Optimizer)\nP_init: initialization matrix (default value is a flat doubly stochastic matrix).\nmax_iter: maximum iterations of the Frank-Wolfe method (default value is 30).\ntol: tolerance for the convergence (default value is 0.1).\n\nOnly for GOAT:\n\nregularization: penalty coefficient in the Sinkhorn algorithm (default value is 100.0).\nmax_iter_sinkhorn: maximum iterations of the Sinkhorn algorithm (default value is 500).\n\nReferences\n\nFAQ: Algorithm 1 of https://arxiv.org/pdf/2111.05366.pdf\nGOAT: Algorithm 3 of https://arxiv.org/pdf/2111.05366.pdf\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.graph_matching_step_size","page":"Algorithms","title":"GraphsOptim.graph_matching_step_size","text":"graph_matching_step_size(A, B, P, Q)\n\nGiven the adjacency matrices A and B, the doubly stochastic matrix P and the direction matrix Q, return the step size of the Frank-Wolfe method for graph matching algorithms.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#Coloring","page":"Algorithms","title":"Coloring","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"fractional_chromatic_number\nfractional_clique_number","category":"page"},{"location":"algorithms/#GraphsOptim.fractional_chromatic_number","page":"Algorithms","title":"GraphsOptim.fractional_chromatic_number","text":"fractional_chromatic_number(g; optimizer)\n\nCompute the fractional chromatic number of a graph.  Gives the same result as fractional_clique_number, though one function may run faster than the other. Beware: this can run very slowly for graphs of any substantial size.\n\nKeyword arguments\n\noptimizer: JuMP-compatible solver (default is HiGHS.Optimizer)\n\nReferences\n\nhttps://mathworld.wolfram.com/FractionalChromaticNumber.html\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.fractional_clique_number","page":"Algorithms","title":"GraphsOptim.fractional_clique_number","text":"fractional_clique_number(g; optimizer)\n\nCompute the fractional clique number of a graph.  Gives the same result as fractional_chromatic_number, though one function may run faster than the other. Beware: this can run very slowly for graphs of any substantial size.\n\nKeyword arguments\n\noptimizer: JuMP-compatible solver (default is HiGHS.Optimizer)\n\nReferences\n\nhttps://mathworld.wolfram.com/FractionalCliqueNumber.html\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#Utils","page":"Algorithms","title":"Utils","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"GraphsOptim.is_binary\nGraphsOptim.is_square\nGraphsOptim.is_stochastic\nGraphsOptim.is_doubly_stochastic\nGraphsOptim.is_permutation_matrix\nGraphsOptim.flat_doubly_stochastic\nGraphsOptim.indvec","category":"page"},{"location":"algorithms/#GraphsOptim.is_binary","page":"Algorithms","title":"GraphsOptim.is_binary","text":"is_binary(x)\n\nCheck if a value is equal to 0 or 1 of its type.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.is_square","page":"Algorithms","title":"GraphsOptim.is_square","text":"is_square(M)\n\nCheck if a matrix is square.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.is_stochastic","page":"Algorithms","title":"GraphsOptim.is_stochastic","text":"is_stochastic(S)\n\nCheck if a matrix is row-stochastic.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.is_doubly_stochastic","page":"Algorithms","title":"GraphsOptim.is_doubly_stochastic","text":"is_doubly_stochastic(D)\n\nCheck if a matrix is row- and column-stochastic.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.is_permutation_matrix","page":"Algorithms","title":"GraphsOptim.is_permutation_matrix","text":"is_permutation_matrix(P)\n\nCheck if a matrix is a permutation matrix.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.flat_doubly_stochastic","page":"Algorithms","title":"GraphsOptim.flat_doubly_stochastic","text":"flat_doubly_stochastic(n)\n\nReturn the barycenter of doubly stochastic matrices J = 𝟏 * 𝟏ᵀ / n.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.indvec","page":"Algorithms","title":"GraphsOptim.indvec","text":"indvec(s, n)\n\nReturn a vector of length n with ones at indices specified by s.\n\n\n\n\n\n","category":"function"}]
}
