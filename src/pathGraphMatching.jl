
"""
pathAlgorithm(G::Matrix{Float64}, H::Matrix{Float64}, ϵ_λ_f::Float64=0.1, ϵ_λ_p::Float64=0.1; kwargs...)

Solve the graph matching problem (or Quadratic Assignment Problem / QAP) using a path-following algorithm.

The algorithm tracks a convex combination (via parameter λ ∈ [0, 1]) between an easily solvable convex relaxation and a concave function that is as hard as the problem. If input matrices differ in size, the smaller matrix is automatically padded with zero rows and columns.

# Arguments
- `G::Matrix{Float64}`: Adjacency or cost matrix of the first graph (n × n).
- `H::Matrix{Float64}`: Adjacency or cost matrix of the second graph (m × m).
- `ϵ_λ_f::Float64=0.1`: Threshold for the maximum normalized change in function value during dynamic step-size control of λ.
- `ϵ_λ_p::Float64=0.1`: Threshold for the normalized permutation matrix change (||P_{new} - P_{opt}|| / √{2n}) during step-size control.

# Keywords
- `dλ_min::Float64=1.0e-5`: Minimum step size for incrementing the path parameter λ.
- `solveQAP::Bool=false`: If `true`, adjusts sign logic to solve a general QAP.
- `return_log::Bool=false`: If `true`, returns a formatted summary log string (runtime, costs, and iterations).
- `return_dataPoints::Bool=false`: If `true`, returns a `NamedTuple` containing trace histories for λ, f_0, f_1, and f_λ.
- `verbose::Bool=false`: Enable detailed console output tracking path-following progress.
- `verbose_FW::Bool=false`: Enable console output for individual Frank-Wolfe optimization steps.

# Returns
- `p_vec::Vector{Int}`: Resulting permutation vector indicating node assignments.
- `log_string::Union{String, Nothing}`: Formatted summary string if `return_log=true`, otherwise `nothing`.
- `dataPoints::Union{NamedTuple, Nothing}`: `NamedTuple` containing `λ_list`, `f0_list`, `f1_list`, and `fλ_list` if `return_dataPoints=true`, otherwise `nothing`.
"""
function pathAlgorithm(
    G::Matrix{Float64},
    H::Matrix{Float64},
    ϵ_λ_f::Float64=0.1,
    ϵ_λ_p::Float64=0.1;
    dλ_min::Float64=1.0e-5,
    solveQAP::Bool=false,
    return_log::Bool=false,
    return_dataPoints::Bool=false,
    verbose::Bool=false,
    verbose_FW::Bool=false,
)
    # extend the smaller matrix by zero rows and columns
    diffSize = size(G, 1) - size(H, 1)

    if diffSize > 0
        # G is larger
        verbose && println(
            "G is larger than H by ", diffSize, " rows and columns. Adding zeros to H."
        )
        H = cat(H, zeros(diffSize, diffSize); dims=(1, 2))
    elseif diffSize < 0
        # H is larger
        verbose && println(
            "H is larger than G by ",
            abs(diffSize),
            " rows and columns. Adding zeros to G.",
        )
        diffSize = abs(diffSize)
        G = cat(G, zeros(diffSize, diffSize); dims=(1, 2))
    end
    m_size = size(G, 1)
    verbose && println("Size of G and H: ", m_size, " x ", m_size)

    t1 = time()

    # allocate fixed space for the gradient matrices so that they don't allocate new space in each calculation
    storage0 = Matrix{Float64}(undef, m_size, m_size)
    storage1 = Matrix{Float64}(undef, m_size, m_size)

    # Start with P as the identity matrix
    p_start = Matrix(1.0I, m_size, m_size)
    lmo = FrankWolfe.BirkhoffPolytopeLMO() #via Hungarian algorithm

    verbose && println("Starting path-following algorithm with λ = 0.0")

    # find initial minimum of F0 ( -F1 for QAP)
    # TODO use Newton instead of FrankWolfe for initialization
    verbose && println("Finding initial minimum of F0 with FrankWolfe")
    if !solveQAP
        init_f = FλForP(0.0, G, H)
        init_∇! = ∇FλForP!(storage0, storage1, 0.0, G, H)
    else
        init_f = FλForP_QAP(0.0, G, H)
        init_∇! = ∇FλForP_QAP!(storage0, storage1, 0.0, G, H)
    end

    p_opt, _ = FrankWolfe.frank_wolfe(
        init_f,
        init_∇!,
        lmo,
        p_start;
        epsilon=1e-8,
        max_iteration=10_000,
        verbose=verbose_FW,
    )

    # change in λ is dynamically adjusted; starts at minimum
    dλ = dλ_min
    # begin with λ=0; iteratively increase up until 1
    λ = 0.0

    # redefine f0, f1 and fλ depending on whether the QAP should be solved or not, s.t. f0 is always convex and f1 is always concave.
    if !solveQAP
        fλNormalizedFinal = fλNormalized
        fλNormalizedFinal = fλ_QAP
    end

    count_iter = 0

    λ_list = Float64[]
    f0_list = Float64[]
    f1_list = Float64[]
    fλ_list = Float64[]
    if return_dataPoints
        push!(λ_list, λ)
        push!(f0_list, f0(p_opt, G, H))
        push!(f1_list, f1(p_opt, G, H))
        push!(fλ_list, fλ(p_opt, λ, G, H))
    end

    verbose && println("λ = ", λ)
    verbose && println()
    while (λ < 1.0)
        count_iter += 1
        # set first possible value for λ_new
        λ_new = λ + dλ

        # calculate local optimum w.r.t. initial λ_new
        verbose && println("   dλ = ", dλ)
        if !solveQAP
            fλ_new_minimize = FλForP(λ_new, G, H)
            ∇fλ_new_minimize = ∇FλForP!(storage0, storage1, λ_new, G, H)
        else
            fλ_new_minimize = FλForP_QAP(λ_new, G, H)
            ∇fλ_new_minimize = ∇FλForP_QAP!(storage0, storage1, λ_new, G, H)
        end
        p_new, _ = frank_wolfe(
            fλ_new_minimize,
            ∇fλ_new_minimize,
            lmo,
            p_opt;
            epsilon=1e-8,
            max_iteration=10_000,
            verbose=verbose_FW,
        )
        p_change_normalized = norm(p_new - p_opt) / sqrt(2 * m_size)

        p_last::Union{Nothing,Matrix{Float64}} = nothing

        # update dλ until criterion is met
        # TODO implemented new stopping criterion. Need to still find out ϵ_f and ϵ_p values from FrankWolfe implementation and calculate ϵ_λ_f and ϵ_λ_p with added input M.
        # d_λ is doubled until one value is larger than it's threshold (or new λ is already 1)
        while abs(
                      fλNormalizedFinal(p_new, λ_new, G, H) -
                      fλNormalizedFinal(p_opt, λ, G, H),
                  ) < ϵ_λ_f &&
                  p_change_normalized < ϵ_λ_p &&
                  λ_new < one(Float64)
            dλ = 2 * dλ
            λ_new = min(λ + dλ, one(Float64))

            verbose && println("   dλ = ", dλ)
            if !solveQAP
                fλ_new_minimize = FλForP(λ_new, G, H)
                ∇fλ_new_minimize = ∇FλForP!(storage0, storage1, λ_new, G, H)
            else
                fλ_new_minimize = FλForP_QAP(λ_new, G, H)
                ∇fλ_new_minimize = ∇FλForP_QAP!(storage0, storage1, λ_new, G, H)
            end
            p_last = p_new
            p_new, _ = frank_wolfe(
                fλ_new_minimize,
                ∇fλ_new_minimize,
                lmo,
                p_opt;
                epsilon=1e-8,
                max_iteration=10_000,
                verbose=verbose_FW,
            )
            p_change_normalized = norm(p_new - p_opt) / sqrt(2 * m_size)
        end

        # if the last while loop's condition is not met (anymore), dλ is one step too large and can be halved once directly
        dλ = max(dλ / 2, dλ_min)
        λ_new = λ + dλ
        verbose && println("   dλ = ", dλ)
        if !isnothing(p_last)
            p_new = p_last
        else
            if !solveQAP
                fλ_new_minimize = FλForP(λ_new, G, H)
                ∇fλ_new_minimize = ∇FλForP!(storage0, storage1, λ_new, G, H)
            else
                fλ_new_minimize = FλForP_QAP(λ_new, G, H)
                ∇fλ_new_minimize = ∇FλForP_QAP!(storage0, storage1, λ_new, G, H)
            end
            p_new, _ = frank_wolfe(
                fλ_new_minimize,
                ∇fλ_new_minimize,
                lmo,
                p_opt;
                epsilon=1e-8,
                max_iteration=10_000,
                verbose=verbose_FW,
            )
        end
        p_change_normalized = norm(p_new - p_opt) / sqrt(2 * m_size)

        # d_λ is halved until both values are smaller than their thresholds (or dλ is already at minimum)
        while (
            abs(fλNormalizedFinal(p_new, λ_new, G, H) - fλNormalizedFinal(p_opt, λ, G, H)) >
            ϵ_λ_f || p_change_normalized > ϵ_λ_p
        ) && dλ > dλ_min
            dλ = max(dλ / 2, dλ_min)
            λ_new = min(λ + dλ, one(Float64))
            verbose && println("   dλ = ", dλ)

            if !solveQAP
                fλ_new_minimize = FλForP(λ_new, G, H)
                ∇fλ_new_minimize = ∇FλForP!(storage0, storage1, λ_new, G, H)
            else
                fλ_new_minimize = FλForP_QAP(λ_new, G, H)
                ∇fλ_new_minimize = ∇FλForP_QAP!(storage0, storage1, λ_new, G, H)
            end
            p_new, _ = frank_wolfe(
                fλ_new_minimize,
                ∇fλ_new_minimize,
                lmo,
                p_opt;
                epsilon=1e-8,
                max_iteration=10_000,
                verbose=verbose_FW,
            )
            p_change_normalized = norm(p_new - p_opt) / sqrt(2 * m_size)
        end
        λ = λ_new
        verbose && println("λ = ", λ)
        verbose && println()
        # criterion is met, λ is set correctly and p_new contans the local optimum w.r.t. the new λ. Set p_opt to p_new for next iteration.

        p_opt = p_new

        if return_dataPoints
            push!(λ_list, λ)
            push!(f0_list, f0(p_opt, G, H))
            push!(f1_list, f1(p_opt, G, H))
            push!(fλ_list, fλ(p_opt, λ, G, H))
        end

        # stop immediately if FrankWolfe arrives at a Permutationmatrix as this is a feasible minimum
        if isPerm(p_opt)
            verbose && println(
                "Found a Permutationmatrix as local optimum, stopping path-following algorithm",
            )
            verbose && println("P:")
            break
        end
    end
    p_vec = permMtV(p_opt)
    verbose && display(p_vec)

    elapsed_time = time() - t1
    verbose && println("Elapsed time: ", elapsed_time, " seconds")

    log_stream = IOBuffer()
    if return_log
        function write_log(msg)
            return println(log_stream, msg) # Schreibt in den Buffer
        end

        write_log("="^60)
        write_log("Results for Graph Matching/QAP")
        write_log("="^60)
        write_log("")
        write_log("ϵ_λ_f: $(ϵ_λ_f)")
        write_log("ϵ_λ_p: $(ϵ_λ_p)")
        write_log("solveQAP: $(solveQAP)")
        write_log("")
        write_log("Runtime: $(elapsed_time) seconds")
        write_log("λ Iterations: $(count_iter)")
        write_log("")
        write_log("Cost:")
        if !solveQAP
            write_log("F0: $(f0(p_opt, G, H))")
            write_log("F1: $(f1(p_opt, G, H))")
        else
            write_log("$(qapVal(p_opt, G, H))")
        end
        write_log("")
        write_log("-"^60)
        write_log("Resulting Matrix P")
        write_log("-"^60)
        write_log(p_vec)
    end
    log_string = return_log ? String(take!(log_stream)) : nothing

    dataPoints = if return_dataPoints
        (; λ_list=λ_list, f0_list=f0_list, f1_list=f1_list, fλ_list=fλ_list)
    else
        nothing
    end

    return p_vec, log_string, dataPoints
end

# returns true if P contains only zeros and ones and false if not
function isPerm(P)
    return all(x -> x == 0.0 || x == 1.0, P)
end

# returns the permutation vector of a permutation matrix P
function permMtV(P)
    return [argmax(row) for row in eachrow(P)]
end

# returns the permutation matrix of a permutation vector P
function permVtM(P)
    return Matrix{Float64}(I(length(P))[P, :])
end

# returns the squared frobenius norm of matrix A
function sqd_frob(A)
    val = norm(A, 2)
    return val^2
end

# returns the diagonal degree matrix of G
# column by column is quicker to go through in julia
function diagonal_degree(G)
    D = zeros(size(G))
    for j in 1:size(D, 1)
        sum = 0.0
        for i in 1:size(D, 1)
            sum += G[i, j]
        end
        D[j, j] = sum
    end
    return D
end

# returns the matrix Δ as stated in the paper
function Δ(G, H)
    D_G = diagonal_degree(G)
    D_H = diagonal_degree(H)

    Δ_G_H = zeros(size(G))
    for i in 1:size(G, 1)
        for j in 1:size(G, 1)
            Δ_G_H[i, j] = D_H[j, j] - D_G[i, i]
        end
    end
    return Δ_G_H .^ 2
end

# returns the laplacian matrix of G
function laplacian(G)
    return diagonal_degree(G) .- G
end

# convex function F0
# algorithm uses only normalized version. This is just for plotting and displaying the correct data.
function f0(P, G, H)
    return sqd_frob(G * P .- P * H)
end

# F0 normalized for values between 0 and 1.
function f0Normalized(P, G, H)
    value = f0(P, G, H)
    return value ./ (sqd_frob(G) + sqd_frob(H))
end

# gradient of F0 but normalized for values between 0 and 1
# save solution value in variable "storage" for space economy
function ∇f0Normalized!(storage, P, G, H)
    value = 2.0 .* ((G^2) * P .- 2.0 .* G * P * H .+ P * (H^2))
    return storage .= value ./ (sqd_frob(G) + sqd_frob(H))
end

# concave function F1.
# algorithm uses only normalized version. This is just for plotting and displaying the correct data.
function f1(P, G, H)
    constantTerm = tr(laplacian(G)^2) + tr(laplacian(H)^2)
    return .-tr(Δ(G, H)' * P) .- 2.0 .* (vec(P)' * vec(laplacian(G) * P * laplacian(H))) +
           constantTerm
end

# F1 normalized for values between 0 and 1.
function f1Normalized(P, G, H)
    value = f1(P, G, H)
    return value ./ (sqd_frob(G) + sqd_frob(H))
end

# gradient of F1 abut normalized for values between 0 and 1
# save solution value in variable "storage" for space economy
function ∇f1Normalized!(storage, P, G, H)
    # the PATH-Algorithm paper has 2.0 in front of the second term, but 4.0 should be correct.
    value = .-Δ(G, H)' .- 4.0 .* laplacian(G) * P * laplacian(H)
    return storage .= value ./ (sqd_frob(G) + sqd_frob(H))
end

# Fλ is convex combination of F0 and F1.
# algorithm uses only normalized version. This is just for plotting and displaying the correct data.
function fλ(P, λ, G, H)
    return (1 - λ) * f0(P, G, H) + λ * f1(P, G, H)
end

# Fλ normalized for values between 0 and 1.
function fλNormalized(P, λ, G, H)
    return (1 - λ) * f0Normalized(P, G, H) + λ * f1Normalized(P, G, H)
end

struct FλForP
    λ::Float64
    G::Matrix{Float64}
    H::Matrix{Float64}
end
# for the FW-algorithm we need a function that takes only P as input.
function (fλ_struct::FλForP)(P)
    return fλNormalized(P, fλ_struct.λ, fλ_struct.G, fλ_struct.H)
end

# function flipped for maximization of the initial function and thus solving QAP
function fλ_QAP(P, λ, G, H)
    return (1 - λ) * (-f1Normalized(P, G, H)) + λ * (-f0Normalized(P, G, H))
end

struct FλForP_QAP
    λ::Float64
    G::Matrix{Float64}
    H::Matrix{Float64}
end
# for the FW-algorithm we need a function that takes only P as input.
function (fλ_struct::FλForP_QAP)(P)
    return fλ_QAP(P, fλ_struct.λ, fλ_struct.G, fλ_struct.H)
end

# gradient of FλNormalized
# save solution value in variable "storage" for space economy
function ∇fλ!(storageλ, storage0, storage1, P, λ, G, H)
    ∇f0Normalized!(storage0, P, G, H)
    ∇f1Normalized!(storage1, P, G, H)
    return storageλ .= (1.0 - λ) .* storage0 .+ λ .* storage1
end
struct ∇FλForP!
    storage0::Matrix{Float64}
    storage1::Matrix{Float64}
    λ::Float64
    G::Matrix{Float64}
    H::Matrix{Float64}
end
# for the FW-algorithm we need a function that takes only P as input.
function (∇fλ_struct::∇FλForP!)(storageλ, P)
    return ∇fλ!(
        storageλ,
        ∇fλ_struct.storage0,
        ∇fλ_struct.storage1,
        P,
        ∇fλ_struct.λ,
        ∇fλ_struct.G,
        ∇fλ_struct.H,
    )
end

# gradient flipped for maximization and solving QAP
# save solution value in variable "storage" for space economy
function ∇fλ_QAP!(storageλ, storage0, storage1, P, λ, G, H)
    ∇f0Normalized!(storage0, P, G, H)
    ∇f1Normalized!(storage1, P, G, H)
    return storageλ .= (1.0 - λ) .* (-storage1) .+ λ .* (-storage0)
end

struct ∇FλForP_QAP!
    storage0::Matrix{Float64}
    storage1::Matrix{Float64}
    λ::Float64
    G::Matrix{Float64}
    H::Matrix{Float64}
end
# for the FW-algorithm we need a function that takes only P as input.
function (∇fλ_struct::∇FλForP_QAP!)(storageλ, P)
    return ∇fλ_QAP!(
        storageλ,
        ∇fλ_struct.storage0,
        ∇fλ_struct.storage1,
        P,
        ∇fλ_struct.λ,
        ∇fλ_struct.G,
        ∇fλ_struct.H,
    )
end

# returns the value of the QAP objective function for a given permutation matrix P and adjacency matrices G and H
function qapVal(P, G, H)
    return tr(G * P * H' * P')
end
