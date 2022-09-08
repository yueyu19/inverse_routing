using JuMP, LinearAlgebra
using Ipopt

"""
Compute entropy-regularized Nash mixed strategies using JuMP and Ipopt
Inputs:
- pa, struct of parameters containing:
    - p, number of players
    - E, incidence matrix
    - n, m = size(E)
    - E_diag, reduced incidence matrix
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

function solve_entropy_routing_jump(pa, b, C, λ, x_init, v_init)
    model = Model(optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-4))
    set_silent(model)

    # variables
    @variable(model, 0 <= x[1:pa.p*pa.m] <= 1)
    @variable(model, v[1:pa.p*(pa.n-1)])

    #todo: warm start x, v
    set_start_value.(x, x_init)
    set_start_value.(v, v_init)

    # slack variables
    @variable(model, p[1:pa.p*pa.m])
    @variable(model, q[1:pa.p*(pa.n-1)])

    # constraints
    for i in 1:pa.p*pa.m
        @NLconstraint(model, x[i] - exp(1/λ * ((pa.E_diag' * v)[i] - b[i] - (C * x)[i]) - ones(pa.p*pa.m,1)[i]) == p[i])
    end

    for j in 1:pa.p*(pa.n-1)
        @NLconstraint(model, pa.s_reduced[j] - (pa.E_diag * x)[j] == q[j])
    end

    # objective
    @objective(model, Min, sum(p[i]*p[i] for i in 1:length(p)) + sum(q[j]*q[j] for j in 1:length(q)))

    optimize!(model)

    println("\t\t\t--- entropy routing solver (JuMP - Ipopt) ---")
    println("\t\t\t termination_status(model) = $(termination_status(model))")
    println("\t\t\t objective_value(model) = $(objective_value(model))")
    println("\t\t\t solve_time(model) = $(solve_time(model))")

    (; x = value.(x), 
       v = value.(v))
end