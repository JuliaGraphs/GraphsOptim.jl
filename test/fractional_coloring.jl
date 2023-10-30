using Graphs
using GraphsOptim
using IterTools
using Test

function kneser_graph(n::Integer, k::Integer)
    ss = collect(subsets(1:n, k))
    return SimpleGraph([isdisjoint(a, b) for a in ss, b in ss])
end

@test fractional_chromatic_number(kneser_graph(8, 3)) ≈ 8 / 3
@test fractional_clique_number(kneser_graph(8, 3)) ≈ 8 / 3
