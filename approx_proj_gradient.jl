include("entropy_routing_solver_jump.jl")
include("routing_solver_jump.jl")
# using IterativeSolvers

# helper function: projection onto D
function proj_D(C, ρ, pa)
    C1 = 0.5 * (C - C')
    for n in 1:pa.p
        C1[pa.m*(n-1)+1:pa.m*n, pa.m*(n-1)+1:pa.m*n] = zeros(pa.m, pa.m)
    end
    
    C2 = 0.5 * (C + C')
    s = eigvals(C2)
    U = eigvecs(C2)
    
    A = real(C1 + U * diagm(vec(max.(s,zeros(length(s))))) * U')
    return ρ / max.(ρ, norm(A)) * A
end


"""
Approximate projected gradient method (measuring function ψ(x)=1/2*||x-x̂||^2)
Inputs:
- p, number of players
- E, incidence matrix
- s, source-sink vector
- x̂, desired Nash solution
- λ, entropy weight
- α, step size
- ϵ, stopping tolerance
- ρ, modification allowance
- max_iter

Returns:
- x: mixed eq strategies for each player on space of links
- b
- C
"""
function approx_proj_grad(p, E, E_diag, s_reduced, x̂, λ, α, ϵ, ρ, max_iter, x_init, b_init, C_init)
    n, m = size(E)
    
    # initiate b, b_plus, C, C_plus
    b = zeros(p*m)
    b_plus = b_init
    C = 0.1*I(p*m)
    C_plus = C_init

    # initiate x, v
    x = x_init
    x_exact = x_init
    v = zeros(p*(n-1))
    v_exact = zeros(p*(n-1))
    
    # data output placeholders
    # ψ_vals = Float64[]
    # violation_metrics = Float64[]
    lambda_vals = Float64[]
    ψ_vals_exact = Float64[]

    println("Starting approx proj grad...")

    for i in 1:max_iter
        # if max(norm(b-b_plus), norm(C-C_plus)) <= ϵ
        if 0.5 * norm(x_exact-x̂)^2 <= ϵ
            println("Converged.")
            break
        elseif i == max_iter
            println("Reached max_iter of $max_iter, break.")
            break
        else
            
            # update b, C
            println("update b, C...\n")
            b = b_plus
            C = C_plus

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

            println("*** Compute (x, v, λ = $λ) => (b, C) ***")
            x, v = solve_entropy_routing(pa(), λ, x, v); println()
            # record current λ value
            push!(lambda_vals, λ)
            
            println("*** Evaluate (b, C) => (x_exact, v_exact) => ψ(x_exact) ***")
            x_exact, v_exact = solve_routing(pa(), x_exact, v_exact)
            @show ψ_x_exact = 0.5 * norm(x_exact-x̂)^2; println()
            # record exact Nash ψ(x)
            push!(ψ_vals_exact, ψ_x_exact)

            # measure ψ(x) val
            # push!(ψ_vals, 0.5 * norm(x-x̂)^2)

            # measure the violation of conditions
            # u = b + C*x - E_diag' * v
            # vio_eqn3 = norm(E_diag*x - s_reduced, Inf)
            # vio_1 = norm(x - max.(0, x-u), Inf)
            # vio_eqn5 = norm(min.(u,0), Inf)
            # vio_eqn2 = norm(u .* x, Inf)
            # vio4 = norm(min.(x, 0), Inf)
            # push!(violation_metrics, max(vio_eqn2, vio_eqn3, vio_eqn5, vio4))
            # push!(violation_metrics, max(vio_eqn3, vio_1))

            # compute D, J
            println("compute D, J...")
            D = diagm(vec(exp.(1/λ * (E_diag'*v-b-C*x) - ones(p*m, 1)))) # dim: 18x18
            J = [I(p*m)+1/λ*D*C     1/λ*D*E_diag'; -E_diag      zeros(p*(n-1), p*(n-1))] # dim: 27x27

            # compute ∇ψ_x
            println("compute ∇ψ_x...")
            ∇ψ_x = x - x̂ # 18x1

            # compute ∇̂ψ_b, ∇̂ψ_C
            println("compute ∇ψ_C...")
            # ∇̂ψ_b = - 1/λ * [D' zeros(p*m, p*(n-1))] * pinv(J)' * [∇ψ_x; zeros(p*(n-1))]
            ∇̂ψ_C = - 1/λ * [D' zeros(p*m, p*(n-1))] * pinv(J)' * [∇ψ_x; zeros(p*(n-1))] * x' # 18x1
            # z = lsmr(J', [∇ψ_x; zeros(p*(n-1))])
            # ∇̂ψ_C = - 1/λ * [D' zeros(p*m, p*(n-1))] * z * x' # 18x1

            # update b_plus, C_plus
            println("compute C_plus...")
            # b_plus = proj_B(b - α * ∇̂ψ_b, 0.1)
            C_plus = proj_D(C - α * ∇̂ψ_C, ρ, pa())
        end
    end

    (;x = x,
      v = v,
      b = b,
      C = C,
    #   ψ_vals = ψ_vals,
    #   violation_metrics = violation_metrics,
      lambda_vals = lambda_vals,
      ψ_vals_exact = ψ_vals_exact)

end