include("routing_solvers/routing_solver_exact_jump.jl")
include("routing_solvers/routing_solver_entropy_nlsolve.jl")
include("routing_solvers/routing_solver_entropy_jump.jl")

# helper function: projection onto B
function proj_B(b)
    return min.(0.1, max.(0, b))
end

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
Approximate projected gradient method
Inputs:
- pa, routing game instance
    - p, number of players
    - E_diag, incidence matrix
    - s_reduced, source-sink vector
    - x̂, desired Nash solution
- λ, entropy weight
- α, step size
- ϵ, stopping tolerance
- ρ, modification allowance
- max_iter
- x_init
- v_init
- b_init
- C_init

Returns:
- x: mixed eq strategies for each player on space of links
- v: value vector/function associated to all nodes
- b: b_i is the nominal price of each link for all i ∈ [p] players
- C: cost adjustment introduced by other players
"""
function approx_proj_grad(pa, λ, α, ϵ, ρ, max_iter, x_init, v_init, b_init, C_init)
    
    # initiate b, b_plus, C, C_plus
    b = zeros(pa.p*pa.m)
    b_plus = b_init
    C = 0.1*I(pa.p*pa.m)
    C_plus = C_init

    # initiate x, v, x_exact, v_exact
    x = x_init
    x_exact = x_init
    v = v_init
    v_exact = v_init

    # initiate ψ_x_exact to check convergence
    ψ_x_exact = Inf

    # data output placeholders
    lambda_vals = Float64[]
    ψ_vals_exact = Float64[]

    for i in 1:max_iter+1
        if ψ_x_exact <= ϵ
            println("\t\tConverged.\n")
            break
        elseif i == max_iter+1
            println("\t\tReached max_iter of $max_iter, break.\n")
            break
        else
            println("\t\t### ITER $i ###\n")

            println("\t\t\t*** Evaluate (b, C) => (x_exact, v_exact) => ψ(x_exact) for iter $i ***")
            x_exact, v_exact = solve_routing(pa, b, C, x_exact, v_exact)
            ψ_x_exact = 0.5 * norm(x_exact-pa.x̂)^2
            println("\t\t\t ψ_x_exact = $ψ_x_exact\n")
            push!(ψ_vals_exact, ψ_x_exact)

            println("\t\t\t*** Compute (x, v, λ = $λ) => (b, C) in iter $i ***")
            # update b, C
            println("\t\t\tupdate b, C...")
            b = b_plus
            C = C_plus

            # call forward solver -- NLsolve or Jump-Ipopt
            if args["forward_solver"] == "nlsolve"
                x, v = solve_entropy_routing_nlsolve(pa, b, C, λ, x, v)
            elseif args["forward_solver"] == "jump"
                x, v = solve_entropy_routing_jump(pa, b, C, λ, x, v)
            end
            push!(lambda_vals, λ)

            # compute D, J
            println("\t\t\tcompute D, J...")
            D = diagm(vec(exp.(1/λ * (pa.E_diag'*v-b-C*x) - ones(pa.p*pa.m, 1)))) # dim: 18x18
            J = [I(pa.p*pa.m)+1/λ*D*C     1/λ*D*pa.E_diag'; -pa.E_diag      zeros(pa.p*(pa.n-1), pa.p*(pa.n-1))] # dim: 27x27

            # compute ∇ψ_x
            println("\t\t\tcompute ∇ψ_x...")
            ∇ψ_x = x - pa.x̂ # 18x1

            # compute ∇̂ψ_b, ∇̂ψ_C
            println("\t\t\tcompute ∇̂ψ_b, ∇ψ_C...")
            ∇̂ψ_b = - 1/λ * [D' zeros(pa.p*pa.m, pa.p*(pa.n-1))] * pinv(J)' * [∇ψ_x; zeros(pa.p*(pa.n-1))]
            ∇̂ψ_C = - 1/λ * [D' zeros(pa.p*pa.m, pa.p*(pa.n-1))] * pinv(J)' * [∇ψ_x; zeros(pa.p*(pa.n-1))] * x' # 18x1

            # update b_plus, C_plus
            println("\t\t\tcompute b_plus, C_plus...\n")
            b_plus = proj_B(b - α * ∇̂ψ_b)
            C_plus = proj_D(C - α * ∇̂ψ_C, ρ, pa)
        end
    end

    (;x = x,
      v = v,
      b = b,
      C = C,
      lambda_vals = lambda_vals,
      ψ_vals_exact = ψ_vals_exact)

end