# flat doubly stochastic matrix, starting point
@doc raw"""
    _flat_doubly_stochastic(n)

Given the dimension ``n``, return the flat doubly stochastic matrix ``J=\frac{1}{n}\mathbb{1}\mathbb{1}^T``.
The flat doubly stochastic matrix is the barycenter of doubly stochastic matrices.
"""
_flat_doubly_stochastic(n::Integer) = ones(n) * ones(n)' / n

# gradient of the objective function
@doc raw"""
    _gradient(A,B,P)

Given the adjacency matrices ``A`` and ``B``, return the gradient evaluated in ``P`` of the function ``f(P)=-\operatorname{trace}(APB^TP^T)``.
"""
function _gradient(A::AbstractMatrix, B::AbstractMatrix, P::AbstractMatrix)
    return -A * P * B' - A' * P * B
end

# norm
"""
    _distance(A,B,P)

Given the adjacency matrices ``A`` and ``B`` and the permutation matrix ``P``, compute the distance between the permuted graphs.
"""
function _distance(A::AbstractMatrix, B::AbstractMatrix, P::AbstractMatrix)
    return norm(A * P - P * B)
end

# solve assignment problem
"""
    _solve_assignment_problem(data; optimizer)

Given a data matrix and a JuMP optimizer, return the integer solution of the assignment problem.
"""
function _solve_assignment_problem(data::AbstractMatrix{U}; optimizer) where {U<:Real}
    M, N = size(data)
    model = Model(optimizer)
    set_silent(model)
    @variable(model, x[m in 1:M, n in 1:N] >= zero(U), integer = true)
    @objective(model, Max, sum(data[m, n] * x[m, n] for n in 1:N, m in 1:M),)
    @constraint(model, [m in 1:M], sum(x[m, :]) == one(U))
    @constraint(model, [n in 1:N], sum(x[:, n]) == one(U))
    optimize!(model)
    P = zeros(Int, M, N)
    P .= value.(x)
    return P
end

"""
    _step_size(A,B,P,Q)

Given ``A,B,P`` and the direction matrix ``Q``, return the step size of the gradient descent method.
"""
function _step_size(
    A::AbstractMatrix, B::AbstractMatrix, P::AbstractMatrix, Q::AbstractMatrix
)
    R = P - Q
    a = -tr(A * R * B' * R')
    b = -tr(A * Q * B' * R' + A * R * B' * Q')
    if a > 0.0
        α = min(max(b / (2 * a), 0.0), 1.0)
    else
        if a + b > 0.0
            α = 0.0
        else
            α = 1.0
        end
    end
    return α
end

"""
    _update_P(P,Q,α)

Given the previous matrix, the direction and the step size, return the updated matrix.
"""
function _update_P(P::AbstractMatrix, Q::AbstractMatrix, α::Real)
    P_new = α * P + (1 - α) * Q
    return P_new
end

"""
    _check_convergence(gradient,P_i,P_new;tol)

Given the gradient and two consecutive step matrices, compute if the method converged.
"""
function _check_convergence(gradient::AbstractMatrix, P_i::AbstractMatrix, P_new::AbstractMatrix; tol::Real)
    if norm(gradient) * (norm(P_i - P_new)) < tol
        return true
    else
        return false
    end
end

"""
    faq(A,B;optimizer,max_iter=30,tol=0.1,init=_flat_doubly_stochastic(size(A)[1]))

Given the adjacency matrix of two graphs, compute the alignment between them. The outputs of the algorithm are:
    - the permutation matrix ``P`` that aligns the two graphs,
    - the distance between the permuted graphs,
    - a boolean that indicates if the algorithm converged.

A JuMP-compatible solver must be provided with the `optimizer` argument.

Ref: [Algorithm 1](https://arxiv.org/pdf/2111.05366.pdf) 

Optional arguments:
    - `max_iter`: maximum of iteration of the gradient descent method, default value is 30.
    - `tol`: tolerance for the convergence, default value is 0.1.
    - `init`: initialization matrix of the gradient descent method, default value is a doubly stochastic matrix.

Example:
```julia-repl
julia> A = rand(10,10);

julia> B = rand(10,10);

julia> faq(a,b;optimizer=HiGHS.Optimizer)
([0 0 … 0 0; 0 0 … 1 0; … ; 0 0 … 0 0; 0 0 … 0 0], 3.9591892429243543, true)
```

"""
function faq(
    A::AbstractMatrix,
    B::AbstractMatrix;
    optimizer,
    max_iter::Integer=30,
    tol::Real=0.1,
    init::AbstractMatrix=_flat_doubly_stochastic(size(A)[1]),
)
    P = init
    converged = false
    for _ in 1:max_iter
        Q = _solve_assignment_problem(_gradient(A, B, P); optimizer=optimizer)
        α_i = _step_size(A, B, P, Q)
        P_new = _update_P(P, Q, α_i)
        converged = _check_convergence(_gradient(A, B, P), P, P_new; tol=tol)
        converged && break
        P = P_new
    end
    P = _solve_assignment_problem(P; optimizer=optimizer)
    return P, _distance(A, B, P), converged
end