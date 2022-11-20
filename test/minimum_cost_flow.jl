using Graphs
using GraphsOptim
using JuMP
using HiGHS
using SparseArrays
using Test

g = DiGraph(Graphs.grid((4, 6)))

vertex_demand = rand(nv(g));
vertex_demand .-= sum(vertex_demand) ./ nv(g);
edge_cost = sparse(src.(edges(g)), dst.(edges(g)), rand(ne(g)));

flow_matrix = minimum_cost_flow(g, vertex_demand, edge_cost)

inflow = vec(sum(flow_matrix; dims=1))
outflow = vec(sum(flow_matrix; dims=2))
@test inflow ≈ outflow + vertex_demand
