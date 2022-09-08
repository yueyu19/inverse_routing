using JuMP, LinearAlgebra
using Ipopt

"""
Compute Exact Nash mixed strategies using JuMP and Ipopt
Inputs:
- pa, struct of parameters containing:
    - p, number of players
    - E, incidence matrix
    - n, m = size(E)
    - E_diag, reduced incidence matrix
    - s_reduced, source-sink vector
- b, where b_i is the nominal price of each link for all i ∈ [p] players
- C, additional cost adjustment introduced by other players
- x_init
- v_init

Returns:
- x: mixed eq strategies for each player on space of links
- v: value vector/function associated to all nodes
"""

function solve_routing(pa, b, C, x_init, v_init)
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
    @constraint(model, b + C * x - pa.E_diag' * v .== u)

    @objective(model, Min, sum(x[i] * u[i] for i in 1:length(x)))

    optimize!(model)

    println("\t\t\t--- exact routing solver ---")
    println("\t\t\t termination_status(model) = $(termination_status(model))")
    println("\t\t\t objective_value(model) = $(objective_value(model))")
    println("\t\t\t solve_time(model) = $(solve_time(model))")

    (; x = value.(x), 
       v = value.(v))
end