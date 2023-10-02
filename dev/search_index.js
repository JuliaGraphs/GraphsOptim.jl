var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"EditURL = \"https://github.com/JuliaGraphs/GraphsOptim.jl/blob/main/README.md\"","category":"page"},{"location":"#GraphsOptim","page":"Home","title":"GraphsOptim","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Dev) (Image: Build Status) (Image: Coverage) (Image: Code Style: Blue)","category":"page"},{"location":"","page":"Home","title":"Home","text":"A package for graph optimization algorithms that rely on mathematical programming.","category":"page"},{"location":"#Getting-started","page":"Home","title":"Getting started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package is still experimental, which is why it is not yet registered. To install it, you need to use the GitHub URL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg; Pkg.add(url=\"https://github.com/JuliaGraphs/GraphsOptim.jl\")","category":"page"},{"location":"#Roadmap","page":"Home","title":"Roadmap","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package only contains a few algorithms, and we would like to add more. New contributors are always welcome: just pick a problem from our roadmap issue, open a pull request following the guidelines, and we will help you get it merged!","category":"page"},{"location":"algorithms/#Algorithms","page":"Algorithms","title":"Algorithms","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"GraphsOptim","category":"page"},{"location":"algorithms/#GraphsOptim","page":"Algorithms","title":"GraphsOptim","text":"GraphsOptim\n\nA package for graph optimization algorithms that rely on mathematical programming.\n\n\n\n\n\n","category":"module"},{"location":"algorithms/#Flow","page":"Algorithms","title":"Flow","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"min_cost_flow\nGraphsOptim.min_cost_flow!","category":"page"},{"location":"algorithms/#GraphsOptim.min_cost_flow","page":"Algorithms","title":"GraphsOptim.min_cost_flow","text":"min_cost_flow(\n    g, vertex_demand, edge_cost, edge_min_capacity, edge_max_capacity;\n    integer, optimizer\n)\n\nCompute a minimum cost flow over a directed graph.\n\nArguments\n\ng::Graphs.AbstractGraph: a directed graph G = (V, E)\nvertex_demand::AbstractVector: a vector in Rⱽ giving the flow requested by each vertex (should be positive for sinks, negative for sources and zero elsewhere)\nedge_cost::AbstractMatrix: a vector in Rᴱ giving the cost of a unit of flow on each edge\nedge_min_capacity::AbstractMatrix: a vector in Rᴱ giving the minimum flow allowed on each edge\nedge_max_capacity::AbstractMatrix: a vector in Rᴱ giving the maximum flow allowed on each edge\n\nKeyword arguments\n\ninteger::Bool: whether the flow should be integer-valued or real-valued\noptimizer: JuMP-compatible solver (default is HiGHS.Optimizer)\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.min_cost_flow!","page":"Algorithms","title":"GraphsOptim.min_cost_flow!","text":"min_cost_flow!(\n    model,\n    g, vertex_demand, edge_cost, edge_min_capacity, edge_max_capacity;\n    var_name, integer\n)\n\nModify a JuMP model by adding the variable, constraints and objective necessary to compute a minimum cost flow over a directed graph.\n\nThe flow variable will be named var_name, see min_cost_flow for details on the other arguments.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"We denote by:","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"f the edge flow variable\nc the edge cost\na and b the min and max edge capacity\nd the vertex demand","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"The objective function is","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"min_f in mathbbR^E sum_(u v) in E c(u v) f(u v)","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"The edge capacity constraint dictates that for all (u v) in E,","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"a(u v) leq f(u v) leq b(u v)","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"The flow conservation constraint with node demand dictates that for all v in V,","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"f^-(v) = d(v) + f^+(v)","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"where the incoming flow f^-(v) and outgoing flow f^+(v) are defined as","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"f^-(v) = sum_u in N^-(v) f(u v) quad textand quad f^+(v) = sum_w in N^+(v) f(v w)","category":"page"},{"location":"algorithms/#Assignment","page":"Algorithms","title":"Assignment","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"danger: Work in progress\nCome back later!","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"min_cost_assignment\nGraphsOptim.min_cost_assignment!","category":"page"},{"location":"algorithms/#GraphsOptim.min_cost_assignment","page":"Algorithms","title":"GraphsOptim.min_cost_assignment","text":"min_cost_assignment(\n    edge_cost;\n    optimizer\n)\n\nCompute a minimum cost assignment over a bipartite graph.\n\nArguments\n\nedge_cost::AbstractMatrix: a matrix in Rᵁˣⱽ giving the cost of matching u ∈ U to v ∈ V\n\nKeyword arguments\n\ninteger::Bool: whether the flow should be integer-valued or real-valued\noptimizer: JuMP-compatible solver (default is HiGHS.Optimizer)\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.min_cost_assignment!","page":"Algorithms","title":"GraphsOptim.min_cost_assignment!","text":"min_cost_assignment!(\n    model,\n    edge_cost;\n    var_name, integer\n)\n\nModify a JuMP model by adding the variable, constraints and objective necessary to compute a minimum cost assignment over a bipartite graph.\n\nThe assignment variable will be named var_name, see min_cost_assignment for details on the other arguments.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#Graph-matching","page":"Algorithms","title":"Graph matching","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"danger: Work in progress\nCome back later!","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"graph_matching\nGraphsOptim.graph_matching_step_size","category":"page"},{"location":"algorithms/#GraphsOptim.graph_matching","page":"Algorithms","title":"GraphsOptim.graph_matching","text":"graph_matching(\n    algo, A, B;\n    optimizer, P_init, max_iter, tol,\n    max_iter_sinkhorn, regularizer\n)\n\nCompute an approximately optimal alignment between two graphs using one of two variants of Frank-Wolfe (FAQ or GOAT)\n\nThe output is a tuple that contains:\n\nthe permutation matrix P defining the alignment;\nthe distance between the permuted graphs;\na boolean indicating if the algorithm converged.\n\nArguments\n\nalgo: allows dispatch based on the choice of algorithm, either FAQ() or GOAT()\nA: the adjacency matrix of the first graph\nB: the adjacency matrix of the second graph\n\nKeyword arguments\n\nFor both algorithms:\n\noptimizer: JuMP-compatible solver (default is HiGHS.Optimizer)\nP_init: initialization matrix (default value is a flat doubly stochastic matrix).\nmax_iter: maximum iterations of the Frank-Wolfe method (default value is 30).\ntol: tolerance for the convergence (default value is 0.1).\n\nOnly for GOAT:\n\nregularization: penalty coefficient in the Sinkhorn algorithm (default value is 100.0).\nmax_iter_sinkhorn: maximum iterations of the Sinkhorn algorithm (default value is 500).\n\nReferences\n\nFAQ: Algorithm 1 of https://arxiv.org/pdf/2111.05366.pdf\nGOAT: Algorithm 3 of https://arxiv.org/pdf/2111.05366.pdf\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.graph_matching_step_size","page":"Algorithms","title":"GraphsOptim.graph_matching_step_size","text":"graph_matching_step_size(A, B, P, Q)\n\nGiven the adjacency matrices A and B, the doubly stochastic matrix P and the direction matrix Q, return the step size of the Frank-Wolfe method for graph matching algorithms.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#Utils","page":"Algorithms","title":"Utils","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"GraphsOptim.is_binary\nGraphsOptim.is_square\nGraphsOptim.is_stochastic\nGraphsOptim.is_doubly_stochastic\nGraphsOptim.is_permutation_matrix\nGraphsOptim.flat_doubly_stochastic","category":"page"},{"location":"algorithms/#GraphsOptim.is_binary","page":"Algorithms","title":"GraphsOptim.is_binary","text":"is_binary(x)\n\nCheck if a value is equal to 0 or 1 of its type.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.is_square","page":"Algorithms","title":"GraphsOptim.is_square","text":"is_square(M)\n\nCheck if a matrix is square.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.is_stochastic","page":"Algorithms","title":"GraphsOptim.is_stochastic","text":"is_stochastic(S)\n\nCheck if a matrix is row-stochastic.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.is_doubly_stochastic","page":"Algorithms","title":"GraphsOptim.is_doubly_stochastic","text":"is_doubly_stochastic(D)\n\nCheck if a matrix is row- and column-stochastic.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.is_permutation_matrix","page":"Algorithms","title":"GraphsOptim.is_permutation_matrix","text":"is_permutation_matrix(P)\n\nCheck if a matrix is a permutation matrix.\n\n\n\n\n\n","category":"function"},{"location":"algorithms/#GraphsOptim.flat_doubly_stochastic","page":"Algorithms","title":"GraphsOptim.flat_doubly_stochastic","text":"flat_doubly_stochastic(n)\n\nReturn the barycenter of doubly stochastic matrices J = 𝟏 * 𝟏ᵀ / n.\n\n\n\n\n\n","category":"function"}]
}
