using Graphs
using GraphsOptim
using JuMP
using HiGHS
using SparseArrays
using Test

g = DiGraph(Graphs.grid((4, 6)))

node_demand = rand(nv(g));
node_demand .-= sum(node_demand) ./ nv(g);
edge_cost = sparse(src.(edges(g)), dst.(edges(g)), rand(ne(g)));

model = Model()
set_optimizer(model, HiGHS.Optimizer)
set_silent(model)

minimum_cost_flow!(model, g, node_demand, edge_cost; flowvar_name=:myflow)

optimize!(model)

@test termination_status(model) == MOI.OPTIMAL
