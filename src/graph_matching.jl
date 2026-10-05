abstract type GraphMatchingAlgorithm end

struct FAQ <: GraphMatchingAlgorithm end

struct GOAT <: GraphMatchingAlgorithm end

"""
    graph_matching_step_size(fP, fQ, cPQ)

Return the step size `α ∈ [0, 1]` maximizing the quadratic objective `f(αP + (1-α)Q) = α² fP + (1-α)² fQ + α(1-α) cPQ` of graph matching along the Frank-Wolfe segment, where `fP = f(P)`, `fQ = f(Q)` and `cPQ = ⟨∇f(P), Q⟩`.
"""
function graph_matching_step_size(fP::Real, fQ::Real, cPQ::Real)
    a = fP + fQ - cPQ
    b = cPQ - 2fQ
    if a < 0
        # concave parabola: clamp its vertex to the segment
        return clamp(-b / (2a), zero(a), one(a))
    else
        # convex or linear: the best endpoint wins
        return fP > fQ ? one(a) : zero(a)
    end
end

"""
    graph_matching_step_size(A, B, P, Q)

Given the adjacency matrices `A` and `B`, the doubly stochastic matrix `P` and the direction matrix `Q`, return the step size of the Frank-Wolfe method for graph matching algorithms, i.e. the `α ∈ [0, 1]` maximizing `tr(A X Bᵀ Xᵀ)` for `X = αP + (1-α)Q`.
"""
function graph_matching_step_size(
    A::AbstractMatrix, B::AbstractMatrix, P::AbstractMatrix, Q::AbstractMatrix
)
    ∇fP = A * P * B' + A' * P * B
    fP = dot(∇fP, P) / 2
    fQ = dot(A * Q * B', Q)
    cPQ = dot(∇fP, Q)
    return graph_matching_step_size(fP, fQ, cPQ)
end

_dense_or_sparse(M::SparseMatrixCSC, ::Type{T}) where {T} = SparseMatrixCSC{T}(M)
_dense_or_sparse(M::AbstractMatrix, ::Type{T}) where {T} = convert(Matrix{T}, M)

"""
    _permute_rows!(X, M, σ)

Fill the dense matrix `X` with `M[σ, :]`.
"""
function _permute_rows!(X::AbstractMatrix, M::AbstractMatrix, σ::AbstractVector{<:Integer})
    @inbounds for j in axes(M, 2), i in axes(M, 1)
        X[i, j] = M[σ[i], j]
    end
    return X
end

function _permute_rows!(X::AbstractMatrix, M::SparseMatrixCSC, σ::AbstractVector{<:Integer})
    fill!(X, zero(eltype(X)))
    σinv = invperm(σ)
    rows, vals = rowvals(M), nonzeros(M)
    @inbounds for j in axes(M, 2), k in nzrange(M, j)
        X[σinv[rows[k]], j] = vals[k]
    end
    return X
end

"""
    _gradient!(G, X, A, B, Bt, P, symmetric)

Overwrite `G` with the gradient `A P Bᵀ + Aᵀ P B` of `P ↦ tr(A P Bᵀ Pᵀ)`, using `X` as a buffer.
"""
function _gradient!(G, X, A, B, Bt, P::AbstractMatrix, symmetric::Bool)
    mul!(X, P, Bt)
    mul!(G, A, X)
    if symmetric
        G .*= 2
    else
        mul!(X, P, B)
        mul!(G, transpose(A), X, true, true)
    end
    return G
end

# For a permutation `Q` with `Q[i, σ[i]] = 1`, we have `Q Bᵀ = Bᵀ[σ, :]` and `Q B = B[σ, :]`
function _gradient!(G, X, A, B, Bt, σ::AbstractVector{<:Integer}, symmetric::Bool)
    _permute_rows!(X, Bt, σ)
    mul!(G, A, X)
    if symmetric
        G .*= 2
    else
        _permute_rows!(X, B, σ)
        mul!(G, transpose(A), X, true, true)
    end
    return G
end

# Barycenter `J = 𝟏𝟏ᵀ / n`: the gradient has rank two and needs no matrix product
function _flat_gradient!(G, A, B)
    n = size(A, 1)
    rA, rB = Vector(vec(sum(A; dims=2))), Vector(vec(sum(B; dims=2)))
    cA, cB = Vector(vec(sum(A; dims=1))), Vector(vec(sum(B; dims=1)))
    G .= (rA .* rB' .+ cA .* cB') ./ n
    return G
end

"""
    _assignment_model(n, optimizer)

Build once the JuMP model of the linear assignment problem (a linear program over the doubly stochastic matrices), whose objective is then updated at every Frank-Wolfe iteration.
"""
function _assignment_model(n::Integer, optimizer)
    model = Model(optimizer)
    set_silent(model)
    min_cost_assignment!(model, zeros(n, n); var_name=:assignment, integer=false)
    return model
end

"""
    _argmax_permutation(score)

Return the row-wise maxima `σ` of `score` if they form a permutation, and `nothing` otherwise.
In the first case, `σ` solves the assignment problem `max ∑ᵢ score[i, σ[i]]`, since the sum of the row-wise maxima is an upper bound.
"""
function _argmax_permutation(score::AbstractMatrix)
    σ = [argmax(view(score, i, :)) for i in axes(score, 1)]
    return isperm(σ) ? σ : nothing
end

"""
    _max_assignment!(model, score)

Return an optimal solution of `max ⟨score, X⟩` over the doubly stochastic matrices `X`.

If `model` is `nothing`, the problem is solved with [`GraphsOptim.linear_sum_assignment`](@ref).
Otherwise, the objective of the JuMP assignment `model` is updated before solving it again, which allows the optimizer to warm start from the previous solution.
"""
function _max_assignment!(model::Union{Nothing,Model}, score::AbstractMatrix)
    σ = _argmax_permutation(score)
    if isnothing(σ) && isnothing(model)
        σ = linear_sum_assignment(score; maximize=true)
    end
    if !isnothing(σ)
        Q = zeros(eltype(score), size(score))
        for i in axes(score, 1)
            Q[i, σ[i]] = 1
        end
        return Q
    end
    x = vec(model[:assignment])
    # Replacing the whole objective at once is much faster than modifying each coefficient
    f = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(-vec(score), index.(x)), 0.0)
    MOI.set(backend(model), MOI.ObjectiveFunction{typeof(f)}(), f)
    optimize!(model)
    @assert termination_status(model) == OPTIMAL
    return reshape(value.(x), size(score))
end

"""
    _permutation_vector(Q)

Return `σ` such that `Q[i, σ[i]] = 1` if `Q` is (numerically) a permutation matrix, and `nothing` otherwise.
"""
function _permutation_vector(Q::AbstractMatrix)
    σ = [argmax(view(Q, i, :)) for i in axes(Q, 1)]
    isperm(σ) || return nothing
    for i in axes(Q, 1)
        isapprox(Q[i, σ[i]], 1; atol=1e-6) || return nothing
    end
    return σ
end

"""
    _distance(P, Q)

Frobenius distance `‖P - Q‖` between two matrices, without allocating `P - Q`.
"""
function _distance(P::AbstractMatrix, Q::AbstractMatrix)
    d = zero(promote_type(eltype(P), eltype(Q)))
    for k in eachindex(P, Q)
        d += abs2(P[k] - Q[k])
    end
    return sqrt(d)
end

function _step_cache(::FAQ, n, optimizer)
    return isnothing(optimizer) ? nothing : _assignment_model(n, optimizer)
end
_step_cache(::GOAT, n, optimizer) = nothing

"""
    _frank_wolfe_step!(algo, P, G, GQ, X, A, B, Bt, symmetric, cache; kwargs...)

Perform one Frank-Wolfe step: compute the search direction `Q` from the gradient `G = ∇f(P)`, store `∇f(Q)` in `GQ`, the optimal step size `α`, and overwrite `P` with `αP + (1-α)Q`.
Return `α` and `‖P - Q‖`.
"""
function _frank_wolfe_step!(::FAQ, P, G, GQ, X, A, B, Bt, symmetric, model; kwargs...)
    n = size(P, 1)
    # Search direction: a vertex Q of the Birkhoff polytope maximizing ⟨∇f(P), Q⟩
    Q = _max_assignment!(model, G)
    σ = _permutation_vector(Q)
    if isnothing(σ)
        _gradient!(GQ, X, A, B, Bt, Q, symmetric)
    else
        # cheaper gradient: Q B and Q Bᵀ are row permutations
        _gradient!(GQ, X, A, B, Bt, σ, symmetric)
    end
    fP = dot(G, P) / 2
    fQ = dot(GQ, Q) / 2
    cPQ = dot(G, Q)
    dist_PQ = _distance(P, Q)
    α = graph_matching_step_size(fP, fQ, cPQ)
    P .= α .* P .+ (1 - α) .* Q
    return α, dist_PQ
end

function _frank_wolfe_step!(
    ::GOAT,
    P,
    G,
    GQ,
    X,
    A,
    B,
    Bt,
    symmetric,
    cache;
    regularization,
    max_iter_sinkhorn,
    tol_sinkhorn,
    kwargs...,
)
    n = size(P, 1)
    # Search direction: the entropic transport plan maximizing ⟨∇f(P), Q⟩,
    # computed with the cost `max(G) - G ≥ 0` (shifting the cost does not change the plan)
    m = maximum(abs, G)
    if m > 0
        Q = sinkhorn(
            ones(eltype(P), n),
            ones(eltype(P), n),
            maximum(G) .- G,
            m / regularization;
            atol=tol_sinkhorn,
            maxiter=max_iter_sinkhorn,
        )
    else
        Q = fill(one(eltype(P)) / n, n, n)
    end
    _gradient!(GQ, X, A, B, Bt, Q, symmetric)
    fP = dot(G, P) / 2
    fQ = dot(GQ, Q) / 2
    cPQ = dot(G, Q)
    dist_PQ = _distance(P, Q)
    α = graph_matching_step_size(fP, fQ, cPQ)
    P .= α .* P .+ (1 - α) .* Q
    return α, dist_PQ
end

"""
    graph_matching(
        algo, A, B;
        optimizer, P_init, max_iter, tol,
        max_iter_sinkhorn, tol_sinkhorn, regularization
    )

Compute an approximately optimal alignment between two graphs using one of two variants of Frank-Wolfe (FAQ or GOAT), by maximizing `tr(A P Bᵀ Pᵀ)` over the doubly stochastic matrices `P`.

The output is a tuple that contains:
- the permutation matrix `P` defining the alignment;
- the distance `‖A P - P B‖` between the permuted graphs;
- a boolean indicating if the algorithm converged.

# Arguments

- `algo`: allows dispatch based on the choice of algorithm, either `FAQ()` or `GOAT()`
- `A`: the adjacency matrix of the first graph
- `B`: the adjacency matrix of the second graph

# Keyword arguments

For both algorithms:

- `optimizer`: JuMP-compatible solver for the linear assignment problems (default is `nothing`, which uses the built-in [`GraphsOptim.linear_sum_assignment`](@ref), much faster for this problem).
- `P_init`: initialization matrix (default value is the flat doubly stochastic matrix).
- `max_iter`: maximum iterations of the Frank-Wolfe method (default value is 30).
- `tol`: the algorithm stops when `‖Pₖ₊₁ - Pₖ‖ / √n < tol` (default value is 0.01).

Only for GOAT:

- `regularization`: penalty coefficient in the Sinkhorn algorithm (default value is 100.0).
- `max_iter_sinkhorn`: maximum iterations of the Sinkhorn algorithm (default value is 1000).
- `tol_sinkhorn`: tolerance on the marginals in the Sinkhorn algorithm in the 1-norm (default value is 0.05).

# References

- FAQ: Algorithm 1 of <https://arxiv.org/pdf/2111.05366.pdf>
- GOAT: Algorithm 3 of <https://arxiv.org/pdf/2111.05366.pdf>
"""
function graph_matching(
    algo::GraphMatchingAlgorithm,
    A::AbstractMatrix,
    B::AbstractMatrix;
    optimizer=nothing,
    P_init::Union{Nothing,AbstractMatrix}=nothing,
    max_iter::Integer=30,
    tol::Real=0.01,
    regularization::Real=100.0,
    max_iter_sinkhorn::Integer=1000,
    tol_sinkhorn::Real=0.05,
)
    n = size(A, 1)
    size(A) == size(B) == (n, n) ||
        throw(DimensionMismatch("A and B must be square matrices of the same size"))
    T = float(promote_type(eltype(A), eltype(B)))
    A = _dense_or_sparse(A, T)
    B = _dense_or_sparse(B, T)
    symmetric = issymmetric(A) && issymmetric(B)
    Bt = symmetric ? B : copy(transpose(B))
    G = Matrix{T}(undef, n, n)  # gradient at P
    GQ = Matrix{T}(undef, n, n)  # gradient at Q
    X = Matrix{T}(undef, n, n)
    cache = _step_cache(algo, n, optimizer)
    # Find a suitable initial position
    if isnothing(P_init)
        P = fill(one(T) / n, n, n)
        _flat_gradient!(G, A, B)
    else
        P = Matrix{T}(P_init)
        _gradient!(G, X, A, B, Bt, P, symmetric)
    end
    # Iterative procedure. The gradient is linear in P, so it can be updated along the
    # iterations with a single gradient evaluation (at the search direction Q) per step.
    converged = false
    for _ in 1:max_iter
        α, dist_PQ = _frank_wolfe_step!(
            algo,
            P,
            G,
            GQ,
            X,
            A,
            B,
            Bt,
            symmetric,
            cache;
            regularization,
            max_iter_sinkhorn,
            tol_sinkhorn,
        )
        G .= α .* G .+ (1 - α) .* GQ
        # Check convergence, knowing that ‖P_new - P‖ = (1 - α) ‖P - Q‖
        converged = (1 - α) * dist_PQ / sqrt(n) < tol
        converged && break
    end
    # Project onto the set of permutation matrices: the assignment LP has integral vertices,
    # and the integer program is only solved if the optimizer returned a non-vertex solution.
    σ = _argmax_permutation(P)
    if isnothing(σ)
        if isnothing(optimizer)
            σ = linear_sum_assignment(P; maximize=true)
        else
            model = cache isa Model ? cache : _assignment_model(n, optimizer)
            σ = _permutation_vector(_max_assignment!(model, P))
            if isnothing(σ)
                P_int = min_cost_assignment(-P; integer=true, optimizer=optimizer)
                σ = _permutation_vector(P_int)
            end
        end
    end
    P = zeros(T, n, n)
    for i in 1:n
        P[i, σ[i]] = one(T)
    end
    # Compute the error ‖A P - P B‖ = ‖A - P B Pᵀ‖
    dist = norm(A - B[σ, σ])
    return P, dist, converged
end
