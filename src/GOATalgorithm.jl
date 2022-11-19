# flat doubly stochastic matrix, starting point
"""
    _flat_doubly_stochastic(n)

Given the dimension ``n``, return the flat doubly stochastic matrix ``J=\frac{1}{n}\mathbb{1}\mathbb{1}^T``.
The flat doubly stochastic matrix is the barycenter of doubly stochastic matrices.
"""
_flat_doubly_stochastic(n::Int64)=ones(n) * ones(n)' / n

# gradient of the objective function
"""
    _gradient(A,B,P_i)

Given the adjacency matrices ``A`` and ``B``, return the gradient evaluated in ``P_i`` of the function ``f(P)=-\operatorname{trace}(APB^TP^T)``.
"""
_gradient(A::AbstractMatrix{U},B::AbstractMatrix{V},P_i::AbstractMatrix{T}) where {U<:Real,V<:Real,T<:Real}  = - A*P_i*B' - A'*P_i*B

# norm
"""
    _distance(A,B,P)

Given the adjacency matrices ``A`` and ``B`` and the permutation matrix ``P``, compute the distance between the permuted graphs.
"""
_distance(A::AbstractMatrix{U},B::AbstractMatrix{V},P::AbstractMatrix{T}) where {U<:Real,V<:Real,T<:Real} = norm(A*P-P*B)

"""
    _solve_transportation_problem(data; optimizer)

Given a data matrix and a JuMP optimizer, return then doubly stochastic solution of the transportation problem.
"""
function _solve_transportation_problem(data::AbstractMatrix{U}; optimizer) where {U<:Real}
    M,N=size(data)
    model = Model(optimizer)
    set_silent(model)
    @variable(model, x[m in 1:M, n in 1:N] >= zero(U),integer=false)
    @objective(
        model,
        Min,
        sum(data[m, n] * x[m, n] for n in 1:N, m in 1:M),
    )
    @constraint(model, [m in 1:M], sum(x[m, :]) == one(U))
    @constraint(model, [n in 1:N], sum(x[:, n]) == one(U))
    optimize!(model)
    Q=zeros(M,N)
    for n in 1:N
        for m in 1:M
            Q[m,n]=value(x[m,n])
        end
    end
    return Q
end

# solve assignment problem
"""
    _solve_assignment_problem(data; optimizer)

Given a data matrix and a JuMP optimizer, return the integer solution of the assignment problem.
"""
function _solve_assignment_problem(data::AbstractMatrix{U}; optimizer) where {U<:Real}
    M,N=size(data)
    model = Model(optimizer)
    set_silent(model)
    @variable(model, x[m in 1:M, n in 1:N] >= zero(U),integer=true)
    @objective(
        model,
        Max,
        sum(data[m, n] * x[m,n] for n in 1:N, m in 1:M),
    )
    @constraint(model, [m in 1:M], sum(x[m, :]) == one(U))
    @constraint(model, [n in 1:N], sum(x[:, n]) == one(U))
    optimize!(model)
    P=zeros(M,N)
    for n in 1:N
        for m in 1:M
            P[m,n]=value(x[m,n])
        end
    end
    return P
end

"""
    _step_size(A,B,P_i,Q_i)

Given ``A,B,P_i`` and the direction matrix ``Q_i``, return the step size of the gradient descent method.
"""
function _step_size(A::AbstractMatrix{U},B::AbstractMatrix{V},P_i::AbstractMatrix{T},Q_i::AbstractMatrix{W}) where {U<:Real,V<:Real,T<:Real,W<:Real}
    R=P_i-Q_i
    a= -tr(A*R*B'*R')
    b= -tr(A*Q_i*B'*R'+A*R*B'*Q_i')
    if a>0.
        α=min(max(-0.5*b/a,0.),1.)
    else
        if a+b>0.
            α=0.
        else
            α=1.
        end
    end
    return α
end

"""
    _update_P(P_i,Q_i,α)

Given the previous matrix, the direction and the step size, return the updated matrix.
"""
function _update_P(P_i::AbstractMatrix{U},Q_i::AbstractMatrix{V},α::Float64) where {U<:Real,V<:Real}
    P_new=α*P_i + (1-α)*Q_i
    return P_new
end

"""
    _check_convergence(gradient,P_i,P_new;tol)

Given the gradient and two consecutive step matrices, compute if the method converged.
"""
function _check_convergence(gradient::AbstractMatrix{U},P_i::AbstractMatrix{V},P_new::AbstractMatrix{T};tol::Float64) where {U<:Real,V<:Real,T<:Real}
    if norm(gradient)*(norm(P_i-P_new))<tol
        return true
    else
        return false
    end
end

"""
    goat(A,B;optimizer,max_iter=30,tol=0.1,init+=_flat_doubly_stochastic(size(A)[1]))

Given the adjacency matrix of two graphs, compute the alignment between them.
Ref: [Algorithm 3](https://arxiv.org/pdf/2111.05366.pdf)

A JuMP-compatible solver must be provided with the `optimizer` argument.

Optional arguments:
    - `max_iter`: maximum of iteration of the gradient descent method, default value is 30.
    - `tol`: tolerance for the convergence, default value is 0.1.
    - `init`: initialization matrix of the gradient descent method, default value is a doubly stochastic matrix.

"""
function goat(A::AbstractMatrix{U},B::AbstractMatrix{T};optimizer,max_iter::Int64=30,tol::Float64=0.1,init::AbstractMatrix{V}=_flat_doubly_stochastic(size(A)[1])) where {U<:Real,V<:Real,T<:Real}
    P_i=init
    converged=false
    for _ in 1:max_iter
        Q_i=_solve_transportation_problem(_gradient(A,B,P_i);optimizer=optimizer)
        α_i=_step_size(A,B,P_i,Q_i)
        P_new=_update_P(P_i,Q_i,α_i)
        converged=_check_convergence(_gradient(A,B,P_i),P_i,P_new;tol=tol)
        converged && break
        P_i=P_new
    end
    P=_solve_assignment_problem(P_i,optimizer=optimizer)
    return P, _distance(A,B,P), converged
end
