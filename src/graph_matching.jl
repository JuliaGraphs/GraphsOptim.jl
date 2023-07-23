abstract type GraphMatchingAlgorithm end

struct FAQ <: GraphMatchingAlgorithm end

struct GOAT <: GraphMatchingAlgorithm end

"""
    graph_matching_step_size(A, B, P, Q)

Given the adjacency matrices `A` and `B`, the doubly stochastic matrix `P` and the direction matrix `Q`, return the step size of the Frank-Wolfe method for graph matching algorithms.
"""
function graph_matching_step_size(
    A::AbstractMatrix, B::AbstractMatrix, P::AbstractMatrix, Q::AbstractMatrix
)
    R = P - Q
    a = -tr(A * R * B' * R')
    b = -tr(A * Q * B' * R' + A * R * B' * Q')
    if a > eps(typeof(a))
        α = min(max(b / (2a), 0.0), 1.0)
    elseif a + b > eps(typeof(a + b))
        α = 0.0
    else
        α = 1.0
    end
    return α
end

function graph_matching_search_direction(::FAQ, ∇f; optimizer, kwargs...)
    return min_cost_assignment(∇f; integer=false, optimizer=optimizer)
end

function graph_matching_search_direction(
    ::GOAT, ∇f::AbstractMatrix{T}; regularization, max_iter_sinkhorn, kwargs...
) where {T}
    m = size(∇f, 1)
    return sinkhorn(
        ones(T, m),
        ones(T, m),
        -∇f,
        -inv(regularization / maximum(abs, ∇f));
        maxiter=max_iter_sinkhorn,
    )
end

"""
    graph_matching(
        algo, A, B;
        optimizer, P_init, max_iter, tol,
        max_iter_sinkhorn, regularizer
    )

Compute an approximately optimal alignment between two graphs using one of two variants of Frank-Wolfe (FAQ or GOAT)

The output is a tuple that contains:
- the permutation matrix `P` defining the alignment;
- the distance between the permuted graphs;
- a boolean indicating if the algorithm converged.

# Arguments

- `algo`: allows dispatch based on the choice of algorithm, either `FAQ()` or `GOAT()`
- `A`: the adjacency matrix of the first graph
- `B`: the adjacency matrix of the second graph

# Keyword arguments

For both algorithms:

- `optimizer`: JuMP-compatible solver (default is `HiGHS.Optimizer`)
- `P_init`: initialization matrix (default value is a flat doubly stochastic matrix).
- `max_iter`: maximum iterations of the Frank-Wolfe method (default value is 30).
- `tol`: tolerance for the convergence (default value is 0.1).

Only for GOAT:

- `regularization`: penalty coefficient in the Sinkhorn algorithm (default value is 100.0).
- `max_iter_sinkhorn`: maximum iterations of the Sinkhorn algorithm (default value is 500).

# References

- FAQ: Algorithm 1 of <https://arxiv.org/pdf/2111.05366.pdf>
- GOAT: Algorithm 3 of <https://arxiv.org/pdf/2111.05366.pdf>
"""
function graph_matching(
    algo::GraphMatchingAlgorithm,
    A::AbstractMatrix,
    B::AbstractMatrix;
    optimizer=HiGHS.Optimizer,
    P_init::AbstractMatrix=flat_doubly_stochastic(size(A, 1)),
    max_iter::Integer=30,
    tol::Real=0.1,
    regularization::Real=100.0,
    max_iter_sinkhorn::Integer=500,
)
    converged = false
    # Find a suitable initial position
    P = P_init
    # Iterative procedure
    for _ in 1:max_iter
        # Compute the gradient of the quadratic objective at the current point 
        ∇f = -A * P * B' - A' * P * B
        # Compute the search direction Q
        Q = graph_matching_search_direction(
            algo, ∇f; optimizer, regularization, max_iter_sinkhorn
        )
        # Compute the step size α
        α = graph_matching_step_size(A, B, P, Q)
        # Update P
        P_new = α * P + (1 - α) * Q
        # Check convergence
        converged = norm(∇f) * (norm(P - P_new)) < tol
        P = P_new
        converged && break
    end
    # Project onto the set of permutation matrices
    P = min_cost_assignment(-P; integer=true, optimizer=optimizer)
    # Compute the error
    dist = norm(A * P - P * B)
    return P, dist, converged
end
