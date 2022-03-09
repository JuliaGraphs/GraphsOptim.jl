"""
    default_weights(g)

Return the binary undirected adjacency matrix of `g`.
"""
function default_weights(g::G) where {G<:AbstractGraph}
    m = spzeros(nv(g), nv(g))
    for e in edges(g)
        m[src(e), dst(e)] = 1
    end
    return m
end
