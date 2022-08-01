using LinearAlgebra
include("entropy_routing_solver.jl")
using Plots

# helper function: projection onto D
function proj_D(C, ρ)
    s = real(eigvals(C))
    U = real(eigvecs(C))
    A = 1/2 * (C - C') + U * diagm(vec(max.(s,zeros(length(s))))) * U'
    return ρ / max(ρ, norm(A)) * A
end

"""
Approximate projected gradient method
Inputs:
- ψ, measuring function, e.g. ψ(x)=1/2*||x-x̂||^2
- λ, entropy weight
- α, step size
- ϵ, stopping tolerance

Returns:
- x: mixed eq strategies for each player on space of links
- b
- C
"""
function approx_proj_grad(p, E, s, λ, α, ϵ)
    n, m = size(E)
    x̂ = [0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0] # 18x1

    # initiate b, b_plus, C, C_plus
    b = zeros(p*m)
    b_plus = 2*ϵ*ones(p*m)
    C = zeros(p*m, p*m)
    C_plus = 2*ϵ*I(p*m)

    ψ_vals = Float64[]

    println("Starting approx proj grad...")

    for i in 1:50
        # update b, C
        b = b_plus
        C = C_plus
    
        # create parameters block and solve for x
        function pa()
            pa_p = p
            pa_E = E
            pa_n, pa_m = n, m
            pa_s = s
            pa_b = b
            pa_C = C
            (; p=pa_p, E=pa_E, n=pa_n, m=pa_m, s=pa_s, b=pa_b, C=pa_C)
        end
        x, v = solve_entropy_routing(pa(), λ)
        # @show 0.5 * norm(x-x̂)^2
        push!(ψ_vals, 0.5 * norm(x-x̂)^2)

        # compute D, J 
        D = diagm(vec(exp.(1/λ * (kron(I(p), E')*v-b-C*x) - ones(p*m, 1)))) # dim: 18x18
        J = [I(p*m)+1/λ*D*C 1/λ*D*kron(I(p), E'); kron(-I(p), E) zeros(p*n, p*n)] # dim: 27x27

        # compute ∇ψ_x
        ∇ψ_x = x - x̂ # 18x1

        # compute ∇̂ψ_b, ∇̂ψ_C
        # ∇̂ψ_b = 1/λ * [D' zeros(p*m, p*n)] * pinv(J)' * [∇ψ_x; zeros(p*n)]
        ∇̂ψ_C = 1/λ * [D' zeros(p*m, p*n)] * pinv(J)' * [∇ψ_x; zeros(p*n)] * x' # 18x1

        # update b_plus, C_plus
        # b_plus = proj_B(b - α * ∇̂ψ_b, 0.1)
        C_plus = proj_D(C - α * ∇̂ψ_C, 0.1)  #ρ=0.1
    end
    println("Finished.")

    (;x = x,
      b = b,
      C = C,
      ψ_vals = ψ_vals)

end

# create routing game instance and parameters
p = 3
E = [1 -1  0  0  1 -1; -1  1  1 -1  0  0; 0  0 -1  1 -1 1]
s = vec([1 0 -1; -1 1 0; 0 -1 1])
λ = 0.01
α = 0.1
ϵ = 0.01

# calling method
x, b, C, ψ_vals = approx_proj_grad(p, E, s, λ, α, ϵ)
plot(ψ_vals)
savefig("ψ_vals_plots.png")