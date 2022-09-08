# include("routing_games.jl")
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
- x_init
- v_init

Returns:
- x: mixed eq strategies for each player on space of links
- v: value vector/function associated to all nodes
"""

function solve_routing(pa, x_init, v_init)
    model = Model(optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-10))
    set_silent(model)

    # variables
    @variable(model, 0 <= x[1:pa.p*pa.m] <= 1)
    @variable(model, v[1:pa.p*(pa.n-1)])
    @variable(model, u[1:pa.p*pa.m] >= 0)

    #todo: warm start x, v
    set_start_value.(x, x_init)
    set_start_value.(v, v_init)

    # constraints
    @constraint(model, pa.E_diag * x .== pa.s_reduced)
    @constraint(model, pa.b + pa.C * x - pa.E_diag' * v .== u)

    @objective(model, Min, sum(x[i] * u[i] for i in 1:length(x)))

    optimize!(model)

    println("--- exact routing solver ---")
    @show termination_status(model)
    @show solve_time(model)

    (; x = value.(x), 
       v = value.(v))
end


# game_name, g, p, E, E_diag, s_reduced, x̂ = grid_graph3_4_players_reduced()
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

# x, v = solve_routing(pa(), 0.5*ones(p*m), zeros(p*(n-1)))

# x, v = solve_routing(pa(), x, v)

# @show x
# @show findall(x->x>0.1, x[1:24])
# @show findall(x->x>0.1, x[25:48])
# @show findall(x->x>0.1, x[49:72])
# @show findall(x->x>0.1, x[73:96])


# @show x - max.(zeros(p*m), x + E_diag'*v - b - C*x)
# @show maximum(x - max.(zeros(p*m), x + E_diag'*v - b - C*x))
# @show minimum(x - max.(zeros(p*m), x + E_diag'*v - b - C*x))

# vio = x - max.(zeros(p*m), x + E_diag'*v - b - C*x)
# @show findall(x->x>0, vio[1:24])
# @show findall(x->x>0, vio[25:48])
# @show findall(x->x>0, vio[49:72])
# @show findall(x->x>0, vio[73:96])