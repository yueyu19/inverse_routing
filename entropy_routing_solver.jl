using NLsolve
using LinearAlgebra

"""
Compute entropy-regularized Nash mixed strategies by Newton's method
Inputs:
- pa, struct of parameters containing:
    - p, number of players
    - E, incidence matrix, where (n, m) = size(E)
    - s, source-sink vector
    - b, 
    - C,
- λ, entropy weight

Returns:
- x: mixed eq strategies for each player on space of links
- v: ?
"""
function solve_entropy_routing(pa, λ)

    function nash!(F, xi)
        # x = xi[1:pa.p*pa.m]
        # v = xi[pa.p*pa.m+1:end]
        # F[1:pa.p*pa.m] = x - exp.(1/λ * (kron(I(pa.p), pa.E') * v - pa.b - pa.C * x) - ones(pa.p*pa.m,1))
        # F[pa.p*pa.m+1:end] = pa.s - kron(I(pa.p), pa.E) * x

        for i in 1:pa.p*pa.m
            F[i] = xi[1:pa.p*pa.m][i] - exp(1/λ * ((kron(I(pa.p), pa.E') * xi[pa.p*pa.m+1:end])[i] - pa.b[i] - (pa.C * xi[1:pa.p*pa.m])[i]) - ones(pa.p*pa.m,1)[i])
        end
        for j in 1:pa.p*pa.n
            F[pa.p*pa.m+j] = pa.s[j] - (kron(I(pa.p), pa.E) * xi[1:pa.p*pa.m])[j]
        end
    end

    x0 = 0.5.*ones(pa.p*pa.m, 1)
    v0 = 0.5.*ones(pa.p*pa.n, 1)
    sol = nlsolve(nash!, [x0; v0])
    # print(sol)

    (;x = sol.zero[1:pa.p*pa.m],
      v = sol.zero[pa.p*pa.m+1:end])
end


"""
creating atomic routing game instance, pa
"""
function pa()
    p = 3
    E = [1 -1  0  0  1 -1; -1  1  1 -1  0  0; 0  0 -1  1 -1 1]
    n, m = size(E)
    s = vec([1 0 -1; -1 1 0; 0 -1 1])
    b = 0.1 * ones(p*m, 1)
    C = zeros(p*m, p*m)
    (; p=p, E=E, n=n, m=m, s=s, b=b, C=C)
end


"""
calling solve_entropy_routing()
"""
x, v = solve_entropy_routing(pa(), 0.01)
@show x