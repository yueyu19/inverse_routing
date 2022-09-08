using LinearAlgebra
using NLsolve
using LineSearches
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

# function solve_entropy_routing(pa, λ)
#     # println("calling solve_entropy_routing...")

#     # R = 1/λ * (pa.p * pa.m * maximum(abs.(pa.C)) + maximum(pa.b))

#     function nash!(F, xi)
#         # x = xi[1:pa.p*pa.m]
#         # v = xi[pa.p*pa.m+1:end]
#         # F[1:pa.p*pa.m] = x - exp.(1/λ * (kron(I(pa.p), pa.E') * v - pa.b - pa.C * x) - ones(pa.p*pa.m,1))
#         # F[pa.p*pa.m+1:end] = pa.s - kron(I(pa.p), pa.E) * x

#         for i in 1:pa.p*pa.m
#             # F[i] = xi[1:pa.p*pa.m][i] - exp(1/λ * ((kron(I(pa.p), pa.E') * xi[pa.p*pa.m+1:end])[i] - pa.b[i] - (pa.C * xi[1:pa.p*pa.m])[i]) - ones(pa.p*pa.m,1)[i])
#             # F[i] = exp(-R) * xi[1:pa.p*pa.m][i] - exp(1/λ * ((kron(I(pa.p), pa.E') * xi[pa.p*pa.m+1:end])[i] - pa.b[i] - (pa.C * xi[1:pa.p*pa.m])[i]) - (R+1) * ones(pa.p*pa.m,1)[i])

#             F[i] = xi[1:pa.p*pa.m][i] - exp(1/λ * ((pa.E_diag' * xi[pa.p*pa.m+1:end])[i] - pa.b[i] - (pa.C * xi[1:pa.p*pa.m])[i]) - ones(pa.p*pa.m,1)[i])

#         end
#         for j in 1:pa.p*(pa.n-1)
#             # F[pa.p*pa.m+j] = pa.s[j] - (kron(I(pa.p), pa.E) * xi[1:pa.p*pa.m])[j]
#             # F[pa.p*pa.m+j] = exp(-R) * (pa.s[j] - (kron(I(pa.p), pa.E) * xi[1:pa.p*pa.m])[j])

#             F[pa.p*pa.m+j] = pa.s_reduced[j] - (pa.E_diag * xi[1:pa.p*pa.m])[j]

#         end
#     end

#     # function j!(H, xi)
#     #     # x = xi[1:pa.p*pa.m]
#     #     # v = xi[pa.p*pa.m+1:end]

#     #     # compute D, J 
#     #     D = diagm(vec(exp.(1/λ * (kron(I(pa.p), pa.E')*xi[pa.p*pa.m+1:end] - pa.b - pa.C*xi[1:pa.p*pa.m]) - ones(pa.p*pa.m, 1)))) # dim: 18x18
        
#     #     J = [I(pa.p*pa.m)+1/λ*D*pa.C 1/λ*D*kron(I(pa.p), pa.E'); kron(-I(pa.p), pa.E) zeros(pa.p*pa.n, pa.p*pa.n)] # dim: 27x27
#     #     for i in 1:pa.p*(pa.m+pa.n), j in 1:pa.p*(pa.m+pa.n)
#     #         H[i, j] = J[i, j]
#     #     end
#     # end

#     x0 = 0.5.*ones(pa.p*pa.m, 1)
#     v0 = 0.5.*ones(pa.p*(pa.n-1), 1)
#     sol = nlsolve(nash!, [x0; v0], autodiff = :forward, show_trace=false, ftol=1e-8, iterations=1000)

#     (;x = sol.zero[1:pa.p*pa.m],
#       v = sol.zero[pa.p*pa.m+1:end])
# end

function solve_entropy_routing(pa, λ, x_init, v_init)
    # println("calling solve_entropy_routing...")

    # R = 1/λ * (pa.p * pa.m * maximum(abs.(pa.C)) + maximum(pa.b))

    function nash!(F, xi)
        # x = xi[1:pa.p*pa.m]
        # v = xi[pa.p*pa.m+1:end]
        # F[1:pa.p*pa.m] = x - exp.(1/λ * (kron(I(pa.p), pa.E') * v - pa.b - pa.C * x) - ones(pa.p*pa.m,1))
        # F[pa.p*pa.m+1:end] = pa.s - kron(I(pa.p), pa.E) * x

        for i in 1:pa.p*pa.m
            # F[i] = xi[1:pa.p*pa.m][i] - exp(1/λ * ((kron(I(pa.p), pa.E') * xi[pa.p*pa.m+1:end])[i] - pa.b[i] - (pa.C * xi[1:pa.p*pa.m])[i]) - ones(pa.p*pa.m,1)[i])
            # F[i] = exp(-R) * xi[1:pa.p*pa.m][i] - exp(1/λ * ((kron(I(pa.p), pa.E') * xi[pa.p*pa.m+1:end])[i] - pa.b[i] - (pa.C * xi[1:pa.p*pa.m])[i]) - (R+1) * ones(pa.p*pa.m,1)[i])

            F[i] = xi[1:pa.p*pa.m][i] - exp(1/λ * ((pa.E_diag' * xi[pa.p*pa.m+1:end])[i] - pa.b[i] - (pa.C * xi[1:pa.p*pa.m])[i]) - ones(pa.p*pa.m,1)[i])

        end
        for j in 1:pa.p*(pa.n-1)
            # F[pa.p*pa.m+j] = pa.s[j] - (kron(I(pa.p), pa.E) * xi[1:pa.p*pa.m])[j]
            # F[pa.p*pa.m+j] = exp(-R) * (pa.s[j] - (kron(I(pa.p), pa.E) * xi[1:pa.p*pa.m])[j])

            F[pa.p*pa.m+j] = pa.s_reduced[j] - (pa.E_diag * xi[1:pa.p*pa.m])[j]

        end
    end

    # function j!(H, xi)
    #     # x = xi[1:pa.p*pa.m]
    #     # v = xi[pa.p*pa.m+1:end]

    #     # compute D, J 
    #     D = diagm(vec(exp.(1/λ * (kron(I(pa.p), pa.E')*xi[pa.p*pa.m+1:end] - pa.b - pa.C*xi[1:pa.p*pa.m]) - ones(pa.p*pa.m, 1)))) # dim: 18x18
        
    #     J = [I(pa.p*pa.m)+1/λ*D*pa.C 1/λ*D*kron(I(pa.p), pa.E'); kron(-I(pa.p), pa.E) zeros(pa.p*pa.n, pa.p*pa.n)] # dim: 27x27
    #     for i in 1:pa.p*(pa.m+pa.n), j in 1:pa.p*(pa.m+pa.n)
    #         H[i, j] = J[i, j]
    #     end
    # end

    # x0 = 0.5.*ones(pa.p*pa.m, 1)
    # v0 = 0.5.*ones(pa.p*(pa.n-1), 1)
    sol = nlsolve(nash!, [x_init; v_init], autodiff = :forward, show_trace=true, ftol=1e-3, iterations=1000, factor=0.1)
    # sol = nlsolve(nash!, [x_init; v_init], autodiff = :forward, show_trace=true, ftol=1e-8, iterations=1000, method=:newton, linesearch=StrongWolfe())

    (;x = sol.zero[1:pa.p*pa.m],
      v = sol.zero[pa.p*pa.m+1:end])
end

game_name, g, p, E, E_diag, s_reduced, x̂ = grid_graph3_4_players_reduced()
n, m = size(E)
b = 0.1*ones(p*m)
C = zeros(p*m, p*m)

x_init = rand(p*m)
v_init = zeros(p*(n-1))

# create parameters block and solve for x
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
x, v = solve_entropy_routing(pa(), 0.01, x_init, v_init)