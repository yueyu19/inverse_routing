include("routing_games.jl")
using JuMP, LinearAlgebra
using Ipopt

"""
Compute entropy-regularized Nash mixed strategies using nonlinear program solver
Inputs:
- pa, struct of parameters containing:
    - p, number of players
    - E, incidence matrix
    - n, m = size(E)
    - E_diag,
    - s_reduced, source-sink vector
- b, where b_i is the nominal price of each link for all i ∈ [p] players
- C, additional cost adjustment introduced by other players
- λ, entropy weight
- x_init
- v_init

Returns:
- x: mixed eq strategies for each player on space of links
- v: value vector/function associated to all nodes
"""

function solve_entropy_routing(pa, λ, x_init, v_init)
    println("create Ipopt model")
    model = Model(optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-10))
    set_silent(model)

    # variables
    println("create variables")
    @variable(model, 0 <= x[1:pa.p*pa.m] <= 1)
    @variable(model, v[1:pa.p*(pa.n-1)])

    #todo: warm start x, v
    # set_start_value.(x, x_init)
    # set_start_value.(v, v_init)

    # slack variables
    println("create slack variables")
    @variable(model, p[1:pa.p*pa.m])
    @variable(model, q[1:pa.p*(pa.n-1)])

    # constraints
    println("set NLconstraints")
    for i in 1:pa.p*pa.m
        @NLconstraint(model, x[i] - exp(1/λ * ((pa.E_diag' * v)[i] - pa.b[i] - (pa.C * x)[i]) - ones(pa.p*pa.m,1)[i]) == p[i])
    end

    for j in 1:pa.p*(pa.n-1)
        @NLconstraint(model, pa.s_reduced[j] - (pa.E_diag * x)[j] == q[j])
    end

    # objective
    println("set objective")
    @objective(model, Min, sum(p[i]*p[i] for i in 1:length(p)) + sum(q[j]*q[j] for j in 1:length(q)))

    println("optimize model...")
    optimize!(model)

    println("--- entropy routing solver ---")
    @show termination_status(model)
    @show objective_value(model)
    @show solve_time(model)

    (; x = value.(x), 
       v = value.(v))
end


# game_name, g, p, E, E_diag, s_reduced, x̂ = grid_graph5_4_players_reduced()
# n, m = size(E)
# b = 0.1*ones(p*m)
# C = zeros(p*m, p*m)

# function pa()
#     pa_p = p
#     pa_E = E
#     pa_E_diag = E_diag
#     pa_n, pa_m = n, m
#     pa_s_reduced = s_reduced
#     pa_b = b
#     pa_C = C
#     (; p=pa_p, E=pa_E, E_diag=pa_E_diag, n=pa_n, m=pa_m, s_reduced=pa_s_reduced, b=pa_b, C=pa_C)
# end

# x, v = solve_entropy_routing(pa(), 0.1, 0.5*ones(p*m), zeros(p*(n-1)))   # λ = 0.003 -> MathOptInterface.LOCALLY_INFEASIBLE
# # @assert findall(x->x>=0.2, x[1:24]) == [9, 13]   # expected [9, 13] for 3x3 4 players
# @show findall(x->x>=0.1, x[1:80])   # expected [33, 37, 41, 45] for 5x5 4 players

# x, v = solve_entropy_routing(pa(), 0.01, x, v)   # λ = 0.003 -> MathOptInterface.LOCALLY_INFEASIBLE
# # @assert findall(x->x>=0.2, x[1:24]) == [9, 13]   # expected [9, 13] for 3x3 4 players
# @show findall(x->x>=0.1, x[1:80])   # expected [33, 37, 41, 45] for 5x5 4 players