"""
    optimal_coloring(g)

Finds a proper coloring using the minimum possible number of colors.  Returns a
Graphs.Coloring object.  Beware: this is an NP-complete problem and so runtime
can in the worst case increase exponentially in the size of the graph.
"""
function optimal_coloring(g::AbstractGraph{T})::Coloring{T} where {T<:Integer}
    # It is not clear whether maximum clique or maximal clique is faster.
    # Maximum clique should generally make coloring faster but of course can
    # itself be time consuming to compute.
    max_clique = argmax(length, maximal_cliques(g))
    #max_clique = independent_set(complement(g), DegreeIndependentSet())
    #@show max_clique

    for χ in length(max_clique):nv(g)
        #println("Trying χ=$χ")
        indexer = reshape(1:(nv(g) * χ), (nv(g), χ))
        cnf = Vector{Int64}[]

        # Seeding the solution on a maximal clique makes it solve faster.
        for (i, v) in enumerate(max_clique)
            push!(cnf, [indexer[v, i]])
        end

        for i in 1:nv(g)
            push!(cnf, indexer[i, :])
        end
        for e in edges(g)
            if src(e) != dst(e)
                for j in 1:χ
                    push!(cnf, [-indexer[src(e), j], -indexer[dst(e), j]])
                end
            end
        end
        sol = solve(cnf)
        if typeof(sol) == Symbol
            if sol == :unsatisfiable
                # continue
            else
                throw(ErrorException("PicoSAT.solve returned unrecognized result status"))
            end
        else
            colormatrix = zeros(Bool, (nv(g), χ))
            # typeassert is needed to get JET.test_package to pass
            colormatrix[filter(i -> i > 0, sol::AbstractArray)] .= true
            colors = map(argmax, eachrow(colormatrix))
            @assert all(colors .>= 1)
            @assert all(colors .<= χ)
            for e in edges(g)
                if src(e) != dst(e)
                    @assert colors[src(e)] != colors[dst(e)]
                end
            end
            return Coloring{T}(χ, colors)
        end
    end
    @assert false
end
