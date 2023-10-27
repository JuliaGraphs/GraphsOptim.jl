"""
    minimum_coloring(
        g, max_nb_colors;
        optimizer
    )

Finds a graph coloring using the smallest possible number of colors, assuming that it will not exceed `max_nb_colors`.

Returns a vector of color indices.

Beware: this is an NP-complete problem and so runtime can in the worst case increase exponentially in the size of the graph.
"""
function minimum_coloring(
    g::AbstractGraph, max_nb_colors::Integer; optimizer=HiGHS.Optimizer
)
    model = Model(optimizer)

    @variable(model, x[1:nv(g), 1:max_nb_colors], Bin)
    @variable(model, y[1:max_nb_colors], Bin)

    @objective(model, Min, sum(y))

    for v in 1:nv(g)
        @constraint(model, sum(x[v, :]) == 1)
    end
    for v in 1:nv(g), c in 1:max_nb_colors
        @constraint(model, x[v, c] <= y[c])
    end
    for c in 1:max_nb_colors
        for e in edges(g)
            u, v = src(e), dst(e)
            @constraint(model, x[u, c] + x[v, c] <= 1)
        end
    end

    set_silent(model)
    optimize!(model)
    # return a vector of color indices
    c = [argmax(value.(x[v, :])) for v in 1:nv(g)]
    return c
end
