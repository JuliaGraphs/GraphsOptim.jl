# Maximum weight clique, maximum independent set and minimum vertex cover: the integer
# programs of GraphsOptim.jl against exact algorithms and heuristics from Graphs.jl and NetworkX.

is_clique(g, c) = all(has_edge(g, u, v) for u in c, v in c if u < v)
is_independent_set(g, s) = !any(has_edge(g, u, v) for u in s, v in s if u < v)
function is_vertex_cover(g, c)
    return (cover = Set(c); all(src(e) in cover || dst(e) in cover for e in edges(g)))
end

pyvertices(s) = pyconvert(Vector{Int}, pylist(s)) .+ 1

"""
    clique_methods(g, w)

Methods computing a maximum weight clique of `g` with vertex weights `w`, all exact.
Each returns a vector of vertices.
"""
function clique_methods(g, w)
    nx = pyimport("networkx")
    h = to_networkx(g; vertex_weights=w)
    return [
        "GraphsOptim.jl" =>
            () -> maximum_weight_clique(g; vertex_weights=w, optimizer=OPTIMIZER),
        "Graphs.jl (all maximal cliques)" =>
            () -> argmax(c -> sum(w[c]), maximal_cliques(g)),
        "NetworkX" => () -> pyvertices(nx.max_weight_clique(h; weight="weight")[0]),
    ]
end

"Methods computing a (maximum or maximal) independent set, as pairs `name => (exact, f)`."
function independent_set_methods(g)
    nx = pyimport("networkx")
    h = to_networkx(g)
    hc = nx.complement(h)
    return [
        "GraphsOptim.jl" =>
            (true, () -> maximum_weight_independent_set(g; optimizer=OPTIMIZER)),
        "NetworkX" => (true, () -> pyvertices(nx.max_weight_clique(hc; weight=nothing)[0])),
        "Graphs.jl (greedy)" => (false, () -> independent_set(g, DegreeIndependentSet())),
        "NetworkX (approximation)" =>
            (false, () -> pyvertices(nx.approximation.maximum_independent_set(h))),
    ]
end

"Methods computing a (minimum or small) vertex cover, as pairs `name => (exact, f)`."
function vertex_cover_methods(g)
    nx = pyimport("networkx")
    h = to_networkx(g)
    return [
        "GraphsOptim.jl" => (true, () -> min_vertex_cover(g; optimizer=OPTIMIZER)),
        "Graphs.jl (greedy)" => (false, () -> vertex_cover(g, DegreeVertexCover())),
        "Graphs.jl (2-approximation)" =>
            (false, () -> vertex_cover(g, RandomVertexCover(); rng=Xoshiro(0))),
        "NetworkX (2-approximation)" =>
            (false, () -> pyvertices(nx.approximation.min_weighted_vertex_cover(h))),
    ]
end

function run_methods!(rows, problem, n, methods, skipped; value, valid)
    for (name, (exact, f)) in methods
        name in skipped && continue
        time, solution = measure(f)
        too_slow(time) && push!(skipped, name)
        ok = !isnothing(solution)
        push!(
            rows,
            (;
                problem,
                n,
                method=name,
                exact,
                time,
                value=ok ? value(solution) : NaN,
                valid=ok && valid(solution),
            ),
        )
    end
end

function run_cliques()
    # compilation
    let g = erdos_renyi(10, 0.5; rng=Xoshiro(0)), w = ones(Int, 10)
        foreach(((_, f),) -> f(), clique_methods(g, w))
        foreach(((_, (_, f)),) -> f(), independent_set_methods(g))
        foreach(((_, (_, f)),) -> f(), vertex_cover_methods(g))
    end

    rows = []
    skipped = Set{String}()
    for n in [20, 40, 60, 80, 100, 150, 200]
        rng = Xoshiro(n)
        g = erdos_renyi(n, 0.5; rng)
        w = rand(rng, 1:10, n)
        methods = [name => (true, f) for (name, f) in clique_methods(g, w)]
        run_methods!(
            rows,
            "maximum weight clique",
            n,
            methods,
            skipped;
            value=c -> sum(w[c]),
            valid=c -> is_clique(g, c),
        )
        @info "Clique" n
    end

    skipped_is, skipped_vc = Set{String}(), Set{String}()
    for n in [25, 50, 100, 200, 400, 800]
        g = erdos_renyi(n, 5 / n; rng=Xoshiro(n))
        run_methods!(
            rows,
            "maximum independent set",
            n,
            independent_set_methods(g),
            skipped_is;
            value=length,
            valid=s -> is_independent_set(g, s),
        )
        run_methods!(
            rows,
            "minimum vertex cover",
            n,
            vertex_cover_methods(g),
            skipped_vc;
            value=length,
            valid=c -> is_vertex_cover(g, c),
        )
        @info "Independent set and vertex cover" n
    end
    write_csv("cliques.csv", rows)
    return nothing
end
