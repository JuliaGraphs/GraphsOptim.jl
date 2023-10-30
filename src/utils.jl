"""
    is_binary(x)

Check if a value is equal to 0 or 1 of its type.
"""
is_binary(x) = x ≈ zero(typeof(x)) || x ≈ one(typeof(x))

"""
    is_square(M)

Check if a matrix is square.
"""
is_square(M::AbstractMatrix) = size(M, 1) == size(M, 2)

"""
    is_stochastic(S)

Check if a matrix is row-stochastic.
"""
function is_stochastic(S::AbstractMatrix)
    for row in eachrow(S)
        if !(sum(row) ≈ one(eltype(S)))
            return false
        end
    end
    return true
end

"""
    is_doubly_stochastic(D)

Check if a matrix is row- and column-stochastic.
"""
is_doubly_stochastic(D::AbstractMatrix) = is_stochastic(D) && is_stochastic(D')

"""
    is_permutation_matrix(P)

Check if a matrix is a permutation matrix.
"""
function is_permutation_matrix(P::AbstractMatrix)
    return is_square(P) && all(is_binary, P) && is_doubly_stochastic(P)
end

"""
    flat_doubly_stochastic(n)

Return the barycenter of doubly stochastic matrices `J = 𝟏 * 𝟏ᵀ / n`.
"""
flat_doubly_stochastic(n::Integer) = ones(n) * ones(n)' / n

"""
    indvec(s, n)

Return a vector of length `n` with ones at indices specified by `s`.
"""
function indvec(s::AbstractVector, n::Integer)
    x = zeros(n)
    x[s] .= 1
    return x
end
