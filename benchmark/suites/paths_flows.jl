# Shortest path, min-cost flow and assignment: the linear programs of GraphsOptim.jl
# against combinatorial algorithms from Graphs.jl, NetworkX and SciPy.

"Random directed graph with average out-degree `d`, plus a Hamiltonian cycle so that it is strongly connected."
function random_digraph(rng::AbstractRNG, n::Integer, d::Real)
    g = erdos_renyi(n, d / n; is_directed=true, rng)
    for i in 1:n
        add_edge!(g, i, mod1(i + 1, n))
    end
    return g
end

"Random integer edge weights in `1:wmax`, as a dense matrix (zero outside of the edges)."
function random_edge_weights(rng::AbstractRNG, g::AbstractGraph, wmax::Integer)
    W = zeros(Int, nv(g), nv(g))
    for e in edges(g)
        W[src(e), dst(e)] = rand(rng, 1:wmax)
    end
    return W
end

path_cost(path, W) = sum(W[path[k], path[k + 1]] for k in 1:(length(path) - 1))

function shortest_path_rows(n)
    rng = Xoshiro(n)
    g = random_digraph(rng, n, 5)
    W = random_edge_weights(rng, g, 100)
    s, t = 1, n ÷ 2
    rows = []
    time, path = measure(() -> shortest_path(g, s, t, W; optimizer=OPTIMIZER))
    push!(
        rows,
        (;
            problem="shortest path",
            n,
            library="GraphsOptim.jl",
            time,
            value=path_cost(path, W),
        ),
    )
    time, path = measure(() -> enumerate_paths(dijkstra_shortest_paths(g, s, W), t))
    push!(
        rows,
        (; problem="shortest path", n, library="Graphs.jl", time, value=path_cost(path, W)),
    )
    return rows
end

function min_cost_flow_rows(n)
    rng = Xoshiro(n)
    g = random_digraph(rng, n, 5)
    cost = random_edge_weights(rng, g, 100)
    capacity = random_edge_weights(rng, g, 20)
    # a tenth of the vertices are sources, another tenth are sinks
    demand = zeros(Int, n)
    vs = randperm(rng, n)
    k = max(1, n ÷ 10)
    for i in 1:k
        supply = rand(rng, 1:10)
        demand[vs[i]] -= supply
        demand[vs[k + i]] += supply
    end
    # the Hamiltonian cycle can carry all the flow, at a high cost, so the problem is feasible
    for i in 1:n
        capacity[i, mod1(i + 1, n)] = sum(abs, demand)
        cost[i, mod1(i + 1, n)] = 100
    end
    rows = []
    time, flow = measure(
        () ->
            min_cost_flow(g, demand, cost, zeros(Int, n, n), capacity; optimizer=OPTIMIZER),
    )
    value = isnothing(flow) ? NaN : sum(cost .* flow)
    push!(rows, (; problem="min-cost flow", n, library="GraphsOptim.jl", time, value))

    nx = pyimport("networkx")
    h = nx.DiGraph()
    for v in 1:n
        h.add_node(v - 1; demand=demand[v])
    end
    for e in edges(g)
        i, j = src(e), dst(e)
        h.add_edge(i - 1, j - 1; weight=cost[i, j], capacity=capacity[i, j])
    end
    time, res = pymeasure(nx.network_simplex, h)
    push!(
        rows,
        (;
            problem="min-cost flow",
            n,
            library="NetworkX",
            time,
            value=pyconvert(Float64, res[0]),
        ),
    )
    return rows
end

function assignment_rows(n; with_jump)
    C = rand(Xoshiro(n), n, n)
    rows = []
    if with_jump
        time, P = measure(() -> min_cost_assignment(C; integer=false, optimizer=OPTIMIZER))
        value = isnothing(P) ? NaN : dot(C, P)
        push!(
            rows, (; problem="assignment", n, library="GraphsOptim.jl (JuMP)", time, value)
        )
    end
    time, σ = measure(() -> GraphsOptim.linear_sum_assignment(C))
    value = sum(C[i, σ[i]] for i in 1:n)
    push!(rows, (; problem="assignment", n, library="GraphsOptim.jl (native)", time, value))
    scipy = pyimport("scipy.optimize")
    np = pyimport("numpy")
    time, res = pymeasure(scipy.linear_sum_assignment, np.ascontiguousarray(C))
    σ = pyconvert(Vector{Int}, res[1]) .+ 1
    push!(
        rows,
        (;
            problem="assignment",
            n,
            library="SciPy",
            time,
            value=sum(C[i, σ[i]] for i in 1:n),
        ),
    )
    return rows
end

function run_paths_flows()
    # compilation
    shortest_path_rows(10)
    min_cost_flow_rows(10)
    assignment_rows(10; with_jump=true)

    rows = []
    for n in [100, 200, 500, 1000, 2000]
        append!(rows, shortest_path_rows(n))
        append!(rows, min_cost_flow_rows(n))
        @info "Paths and flows" n
    end
    jump_ok = true
    for n in [10, 20, 50, 100, 200, 500, 1000]
        new_rows = assignment_rows(n; with_jump=jump_ok)
        # stop using JuMP once it hits the time limit
        jump_ok = jump_ok && !isnan(first(new_rows).time)
        append!(rows, new_rows)
        @info "Assignment" n
    end
    write_csv("paths_flows.csv", rows)
    return nothing
end
