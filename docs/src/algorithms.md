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

We denote:

- ``G = (V, E)`` the graph we work on with vertex set ``V`` and edge set ``E``
- ``s`` the source vertex
- ``t`` the target/destination vertex
- ``c`` the edge cost function, i.e. the cost to travel through a specific edge

The objective function is

```math
\min_{x \in \mathbb{R}^E} \sum_{(u, v) \in E} c(u, v) x(u, v)
```

The variable $x(u,v)$ is an edge selection variable. $x(u,v) = 1$ means that edge $(u, v)$ is selected in the path and $x(u, v)=  0$ otherwise.

The linear program is constrained by

1. 0 to 1 selection for all edges:
   ```math
   0 \leq x(u, v) \leq 1
   ```

2. Path goes out from source $s$
   ```math
   \sum_{(s, v) \in E} x(s, v) - \sum_{(v, s) \in E} x(v, s) = 1
   ```

3. Path ends at target $t$
   ```math
   \sum_{(t, v) \in E} x(t, v) + \sum_{(v, t) \in E} x(v, t) = 1
   ```
   
4. All middle vertices must have equal incoming and outgoing paths, i.e. $\forall v \notin \lbrace s, t\rbrace$
   ```math
   \sum_{(u, v) \in E} x(u, v) - \sum_{v, w} x(v, w) = 0
   ```

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
