"""
    faq(
        A, B;
        max_iter, tol, P_init, optimizer
    )

Compute an approximately optimal alignment between two graphs.

The output is a tuple that contains:
- the permutation matrix `P` defining the alignment;
- the distance between the permuted graphs;
- a boolean indicating if the algorithm converged.

Ref: Algorithm 1 of <https://arxiv.org/pdf/2111.05366.pdf>

# Arguments

- `A`: the adjacency matrix of the first graph
- `B`: the adjacency matrix of the second graph

# Keyword arguments

- `max_iter`: maximum iterations of the gradient descent method (default value is 30).
- `tol`: tolerance for the convergence (default value is 0.1).
- `init`: initialization matrix of the gradient descent method (default value is a flat doubly stochastic matrix).
- `optimizer`: JuMP-compatible solver (default is `HiGHS.Optimizer`)
"""
function faq(
    A::AbstractMatrix,
    B::AbstractMatrix;
    max_iter::Integer=30,
    tol::Real=0.1,
    P_init::AbstractMatrix=flat_doubly_stochastic(size(A, 1)),
    optimizer=HiGHS.Optimizer,
)
    converged = false
    # Find a suitable initial position
    P = P_init
    # Iterative procedure
    for _ in 1:max_iter
        # Compute the gradient of the quadratic objective at the current point 
        ∇f = -A * P * B' - A' * P * B
        # Compute the search direction Q
        Q = min_cost_assignment(∇f; integer=false, optimizer=optimizer)
        # Compute the step size α
        α = _step_size(A, B, P, Q)
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

"""
    _step_size(A, B, P, Q)

Given the adjacency matrices `A` and `B`, the doubly stochastic matrix `P` and the direction matrix `Q`, return the step size of the gradient descent method.
"""
function _step_size(
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

function goat(
    A::AbstractMatrix,
    B::AbstractMatrix;
    max_iter::Integer=30,
    tol::Real=0.1,
    P_init::AbstractMatrix=flat_doubly_stochastic(size(A, 1)),
    regulizer::Real=100.0,
    max_iter_sinkhorn::Integer=500,
    optimizer=HiGHS.Optimizer,
)
    m = size(A, 1)
    converged = false
    # Find a suitable initial position
    P = P_init
    # Iterative procedure
    for _ in 1:max_iter
        ∇f = -A * P * B' - A' * P * B
        Q = sinkhorn(ones(m), ones(m), - ∇f, -1/(regulizer/ maximum(abs.(∇f))); max_iter = max_iter_sinkhorn)
        α = _step_size(A,B,P,Q)
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

   