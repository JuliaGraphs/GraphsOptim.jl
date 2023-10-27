using Graphs
using GraphsOptim
using IterTools
using Test

function queens_graph(n::Integer, m::Integer)
    g = SimpleGraph(n * m)
    for (ix, iy, jx, jy) in product(1:n, 1:m, 1:n, 1:m)
        dx = ix - jx
        dy = iy - jy
        if dx != 0 || dy != 0
            if dx == 0 || dy == 0 || dx == dy || dx == -dy
                add_edge!(g, (ix - 1) * m + iy, (jx - 1) * m + jy)
            end
        end
    end
    return g
end

queens_graph(n::Integer) = queens_graph(n, n)

function kneser_graph(n::Integer, k::Integer)
    ss = collect(subsets(1:n, k))
    return SimpleGraph([isdisjoint(a, b) for a in ss, b in ss])
end

# https://oeis.org/A088202

nb_colors(c) = length(unique(c))

@test minimum_coloring(queens_graph(1), 10) |> nb_colors == 1
@test minimum_coloring(queens_graph(2), 10) |> nb_colors == 4
@test minimum_coloring(queens_graph(3), 10) |> nb_colors == 5
@test minimum_coloring(queens_graph(4), 10) |> nb_colors == 5
@test minimum_coloring(queens_graph(5), 10) |> nb_colors == 5
@test minimum_coloring(queens_graph(6), 10) |> nb_colors == 7

@test minimum_coloring(kneser_graph(11, 4), 10) |> nb_colors == 11 - 2 * 4 + 2
