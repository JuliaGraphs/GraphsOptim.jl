# Algorithms

```@docs
GraphsOptim
```

## Flow

```@docs
min_cost_flow
GraphsOptim.min_cost_flow!
```

We denote by:

- ``f`` the edge flow variable
- ``c`` the edge cost
- ``a`` and ``b`` the min and max edge capacity
- ``d`` the vertex demand

The objective function is

```math
\min_{f \in \mathbb{R}^E} \sum_{(u, v) \in E} c(u, v) f(u, v)
```

The edge capacity constraint dictates that for all $(u, v) \in E$,

```math
a(u, v) \leq f(u, v) \leq b(u, v)
```

The flow conservation constraint with node demand dictates that for all $v \in V$,

```math
f^-(v) = d(v) + f^+(v)
```

where the incoming flow $f^-(v)$ and outgoing flow $f^+(v)$ are defined as

```math
f^-(v) = \sum_{u \in N^-(v)} f(u, v) \quad \text{and} \quad f^+(v) = \sum_{w \in N^+(v)} f(v, w)
```

## Shortest Path

```@docs
shortest_path
GraphsOptim.shortest_path!
```

A special case of minimum cost flow without edge capacities, and where vertex demands are $0$ everywhere except at the source ($-1$) and target ($+1$).

## Assignment

!!! danger "Work in progress"
    Come back later!

```@docs
min_cost_assignment
GraphsOptim.min_cost_assignment!
```

## Minimum Vertex Cover

```@docs
min_vertex_cover
GraphsOptim.min_vertex_cover!
```

Finds a subset $S \subset V$ of vertices of an undirected graph $G = (V,E)$ such that $\forall (u,v) \in E: u \in S \lor v \in S$

## Maximum weight clique

```@docs
maximum_weight_clique
GraphsOptim.maximum_weight_clique!
```

A *clique* is a subset $S \subset V$ of vertices of an undirected graph $G = (V,E)$ such that $\forall (u,v) \in S: (u, v) \in E$. We search for the clique maximizing the total weight of selected vertices.

## Maximum Weight Independent Set

```@docs
maximum_weight_independent_set
GraphsOptim.maximum_weight_independent_set!
```

Finds a subset $S \subset V$ of vertices of maximal weight of an undirected graph $G = (V,E)$ such that $\forall (u,v) \in E: u \notin S \lor v \notin S$.

## Graph matching

!!! danger "Work in progress"
    Come back later!

```@docs
graph_matching
GraphsOptim.graph_matching_step_size
```

## Graph edit distance
```@docs
GraphsOptim.edit_distance
GraphsOptim.edit_distance!
```

### Mathematical Details

The formulations implemented here are described in [`D'ascenzo, Andrea, et al. "Enhancing Graph Edit Distance Computation: Stronger and Orientation-based ILP Formulations." Proceedings of the VLDB Endowment 18.11 (2025): 4737-4749](https://www.vldb.org/pvldb/vol18/p4737-d%27ascenzo.pdf)`. The paper provides a more detailed description of the graph edit distance problem.

In the graph edit distance problem, we search for a minimum cost *edit path* between
two graphs $G$ and $H$, that is a sequence of edit operations transforming $G$ to $H$. Valid edit operations include inserting, deleting, or relabeling vertices and edges.

This formulation of the problem is easy to visualize, but is not well suited to an integer-programming implementation. We thus use an equivalent definition: a *node map*
is a relation $\pi \subset V_{G+ \epsilon} \times V_{H + \epsilon}$ on the vertex sets
augmented by $\epsilon$, in which 

* each node $v \in V_G$ is mapped to exactly one element of $V_{H + \epsilon}$
* each node $w \in V_H$ has exactly one preimage in $V_{G + \epsilon}$

Mapping a vertex of $G$ to $\epsilon$ represents a deletion, whereas mapping $\epsilon$ to a vertex of $H$ represents an insertion. For metric cost functions this is equivalent to the edit path formulation (and any cost function can be transformed into an equivalent metric cost function).

This version of the problem is much better suited for integer programming, and is the
foundation for all formulations implemented here. All formulations follow the same basic principle: they define binary variables that encode vertex and edge mappings, impose constraints that ensure valid mappings, and add *topology constraints* that link the vertex and edge variables according to the input graphs.
The formulations differ in the exact variable layout and more importantly the topology
constraints used. As seen in the paper mentioned above, the different topology constraints
not only provide massive speedups in practice, but  also yield relaxations with different bounds.

The best formulation both in theory and practice orients the graphs to provide stricter
topology constraints. The graph $G$ is oriented in a canonical way, while in $H$ there are
forward and backward edges for each undirected edge in the input. This allows us to require edge mappings to preserve orientation.

### Function documentation for advanced use
```@docs
GraphsOptim.EditCosts
GraphsOptim.get_default_edit_costs
GraphsOptim.validate_cost_function
GraphsOptim.Formulation
GraphsOptim.F1
GraphsOptim.F1prime
GraphsOptim.F1plus
GraphsOptim.F2minus
GraphsOptim.F2
GraphsOptim.F2plus
GraphsOptim.FORI
```

## Coloring

```@docs
fractional_chromatic_number
fractional_clique_number
```

## Utils

```@docs
GraphsOptim.is_binary
GraphsOptim.is_square
GraphsOptim.is_stochastic
GraphsOptim.is_doubly_stochastic
GraphsOptim.is_permutation_matrix
GraphsOptim.flat_doubly_stochastic
GraphsOptim.indvec
```
