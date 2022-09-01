using JuMP, Ipopt
using LinearAlgebra
include("routing_games.jl")

"""
Compute entropy-regularized Nash mixed strategies by Newton's method
Inputs:
- pa, struct of parameters containing:
    - p, number of players
    - E, incidence matrix, where (n, m) = size(E)
    - s, source-sink vector
    - b, where b_i is the nominal price of each link for all i ∈ [p] players
    - C, additional cost adjustment introduced by other players
- λ, entropy weight

Returns:
- x: mixed eq strategies for each player on space of links
- v: value vector/function associated to all nodes
"""

function solve_entropy_routing(pa, λ)
    model = Model(Ipopt.Optimizer)
    set_silent(model)

    @variable(model, 0 <= x[1:pa.p*pa.m] <= 1)
    @variable(model, v[1:pa.p*(pa.n-1)])

    # x - exp.(1/λ * (kron(I(pa.p), pa.E') * v - pa.b - pa.C * x) - ones(pa.p*pa.m,1)) .== 0
    for i in 1:pa.p*pa.m
        @NLconstraint(model, x[i] - exp(1/λ * ((pa.E_diag' * v)[i] - pa.b[i] - (pa.C * x)[i]) - ones(pa.p*pa.m,1)[i]) == 0)
    end

    # pa.s - kron(I(pa.p), pa.E) * x .== 0
    for j in 1:pa.p*(pa.n-1)
        @NLconstraint(model, pa.s_reduced[j] - (pa.E_diag * x)[j] == 0)
    end

    # @NLobjective(model, Min, norm([x - exp.(1/λ * (kron(I(p), E') * v - b - C * x) - ones(p*m,1)); pa.s - kron(I(p), E) * x]))

    # print(model)
    optimize!(model)
    @show termination_status(model)

    (; x = value.(x), 
       v = value.(v))
end


game_name, g, p, E, E_diag, s_reduced, x̂ = grid_graph3_4_players_reduced()
n, m = size(E)
b = 0.1*ones(p*m)
C = zeros(p*m, p*m)

function pa()
    pa_p = p
    pa_E = E
    pa_E_diag = E_diag
    pa_n, pa_m = n, m
    pa_s_reduced = s_reduced
    pa_b = b
    pa_C = C
    (; p=pa_p, E=pa_E, E_diag=pa_E_diag, n=pa_n, m=pa_m, s_reduced=pa_s_reduced, b=pa_b, C=pa_C)
end

x, v = solve_entropy_routing(pa(), 0.004)   # λ = 0.003 -> MathOptInterface.LOCALLY_INFEASIBLE
@show findall(x->x>=0.2, x[1:24])   # expected [9, 13]