# Graph matching: FAQ and GOAT against graspologic.

"""
    correlated_er(rng, n, p, ρ)

Sample a pair of `ρ`-correlated Erdős–Rényi graphs with edge probability `p`, and shuffle the vertices of the second one.
Return the adjacency matrices `A`, `B` and the ground truth `σ` such that vertex `i` of `A` corresponds to vertex `σ[i]` of `B`.
"""
function correlated_er(rng::AbstractRNG, n::Integer, p::Real, ρ::Real)
    A, B = zeros(n, n), zeros(n, n)
    for j in 1:n, i in 1:(j - 1)
        a = rand(rng) < p
        b = a ? rand(rng) < p + ρ * (1 - p) : rand(rng) < p * (1 - ρ)
        A[i, j] = A[j, i] = a
        B[i, j] = B[j, i] = b
    end
    π = randperm(rng, n)
    return A, B[π, π], invperm(π)
end

"Number of edges present in exactly one of the two aligned graphs, i.e. `‖A - P B Pᵀ‖² / 2`."
disagreements(A, B, σ) = round(Int, sum(abs2, A - B[σ, σ]) / 2)

function graph_matching_rows(algo, A, B, truth; seed, min_time=1.0)
    name = string(nameof(typeof(algo)))
    graspologic = pyimport("graspologic.match")
    np = pyimport("numpy")
    t_jl, (P, _, _) = measure(() -> graph_matching(algo, A, B); min_time)
    σ_jl = [argmax(view(P, i, :)) for i in axes(P, 1)]
    t_py, res = pymeasure(
        graspologic.graph_match,
        np.ascontiguousarray(A),
        np.ascontiguousarray(B);
        transport=algo isa GOAT,
        rng=seed,
        min_time,
    )
    σ_py = pyconvert(Vector{Int}, res.indices_B[np.argsort(res.indices_A)]) .+ 1
    return [
        (;
            algorithm=name,
            library,
            time,
            accuracy=mean(σ .== truth),
            disagreements=disagreements(A, B, σ),
            disagreements_truth=disagreements(A, B, truth),
        ) for
        (library, time, σ) in (("GraphsOptim.jl", t_jl, σ_jl), ("graspologic", t_py, σ_py))
    ]
end

function run_graph_matching()
    # compilation
    let (A, B, truth) = correlated_er(Xoshiro(0), 20, 0.3, 0.9)
        graph_matching_rows(FAQ(), A, B, truth; seed=0, min_time=0)
        graph_matching_rows(GOAT(), A, B, truth; seed=0, min_time=0)
    end

    # Running time as the number of vertices grows
    scaling = []
    for n in [50, 100, 200, 500, 1000], algo in (FAQ(), GOAT())
        A, B, truth = correlated_er(Xoshiro(n), n, 0.1, 0.9)
        for row in graph_matching_rows(algo, A, B, truth; seed=n)
            push!(scaling, (; n, row...))
        end
        @info "Graph matching" n algo
    end
    write_csv("graph_matching_scaling.csv", scaling)

    # Matching quality as the correlation between the graphs varies
    quality = []
    for ρ in 0.5:0.1:1.0, seed in 1:10, algo in (FAQ(), GOAT())
        A, B, truth = correlated_er(Xoshiro(seed), 150, 0.1, ρ)
        for row in graph_matching_rows(algo, A, B, truth; seed, min_time=0)
            push!(quality, (; ρ, seed, row...))
        end
    end
    write_csv("graph_matching_quality.csv", quality)
    return nothing
end
