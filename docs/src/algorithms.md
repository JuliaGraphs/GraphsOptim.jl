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

## Minimum Vertex Cover

```@docs
min_vertex_cover
```
Minumum vertex cover algorithm finds a subset $S \subset V$ of vertices of an undirected graph $G = (V,E)$ such that $\forall (u,v) \in E: u \in S \lor v \in S$


## Assignment

!!! danger "Work in progress"
    Come back later!

```@docs
min_cost_assignment
GraphsOptim.min_cost_assignment!
```

## Graph matching

!!! danger "Work in progress"
    Come back later!

```@docs
graph_matching
GraphsOptim.graph_matching_step_size
```

## Utils

```@docs
GraphsOptim.is_binary
GraphsOptim.is_square
GraphsOptim.is_stochastic
GraphsOptim.is_doubly_stochastic
GraphsOptim.is_permutation_matrix
GraphsOptim.flat_doubly_stochastic
```
