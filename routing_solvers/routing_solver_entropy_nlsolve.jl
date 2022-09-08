using LinearAlgebra
using NLsolve

"""
Compute entropy-regularized Nash mixed strategies using NLsolve
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

function solve_entropy_routing_nlsolve(pa, b, C, λ, x_init, v_init)
    println("\t\t\t--- entropy routing solver (NLsolve) ---")

    function nash!(F, xi)
        # x = xi[1:pa.p*pa.m]
        # v = xi[pa.p*pa.m+1:end]
        # F[1:pa.p*pa.m] = x - exp.(1/λ * (kron(I(pa.p), pa.E') * v - b - pa.C * x) - ones(pa.p*pa.m,1))
        # F[pa.p*pa.m+1:end] = pa.s - kron(I(pa.p), pa.E) * x

        for i in 1:pa.p*pa.m
            F[i] = xi[1:pa.p*pa.m][i] - exp(1/λ * ((pa.E_diag' * xi[pa.p*pa.m+1:end])[i] - b[i] - (C * xi[1:pa.p*pa.m])[i]) - ones(pa.p*pa.m,1)[i])
        end
        for j in 1:pa.p*(pa.n-1)
            F[pa.p*pa.m+j] = pa.s_reduced[j] - (pa.E_diag * xi[1:pa.p*pa.m])[j]
        end
    end

    sol = nlsolve(nash!, [x_init; v_init], autodiff = :forward, show_trace=false, ftol=1e-4, iterations=1000)

    (;x = sol.zero[1:pa.p*pa.m],
      v = sol.zero[pa.p*pa.m+1:end])
end