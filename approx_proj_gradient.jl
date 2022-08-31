include("entropy_routing_solver_nlsolve.jl")
using IterativeSolvers

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

function vio(b, C, x, v)
    u = b + C*x - kron(I(p), E')*v
    vio_eqn5 = norm(min.(u,0), Inf)
    vio_eqn2 = norm(u .* x, Inf)
    vio_eqn3 = norm(kron(I(p), E)*x - s, Inf)
    vio4 = norm(min.(x, 0), Inf)
    return max(vio_eqn2, vio_eqn3, vio_eqn5, vio4)
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
function approx_proj_grad(p, E, s, x̂, λ, α, ϵ, ρ, max_iter, x_init, b_init, C_init)
    n, m = size(E)
    
    # initiate b, b_plus, C, C_plus AND x
    b = zeros(p*m)
    b_plus = b_init
    C = 0.1*I(p*m)
    C_plus = C_init
    x = x_init

    v = zeros(p*n)
    
    # data output placeholders
    ψ_vals = Float64[]
    violation_metrics = Float64[]
    lambda_vals = Float64[]
    ∇̂ψ_C_norm = Float64[]
    D_norm = Float64[]
    J_norm = Float64[]
    pinv_J_norm = Float64[]
    F_norm = Float64[]

    println("Starting approx proj grad...")

    for i in 1:max_iter
        # if max(norm(b-b_plus), norm(C-C_plus)) <= ϵ || 0.5 * norm(x-x̂)^2 <= ϵ
        # if max(norm(b-b_plus), norm(C-C_plus)) <= ϵ
        # if 0.5 * norm(x-x̂)^2 <= ϵ
        # if vio(b, C, x, v) <= ϵ
        #     println("Converged.")
        #     break
        # elseif i == max_iter
        #     println("Reached max_iter of $max_iter, break.")
        #     break
        # else
            
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
            
            # record norm of F -- violation of nonlinear equations
            F = [x - exp.(1/λ * (kron(I(p), E') * v - b - C * x) - ones(p*m,1)); pa.s - kron(I(p), E) * x]
            push!(F_norm, norm(F))

            # measure ψ(x) val
            push!(ψ_vals, 0.5 * norm(x-x̂)^2)

            # measure the violation of conditions
            push!(violation_metrics, vio(b, C, x, v))

            # record current λ value
            push!(lambda_vals, λ)
            
            # compute D, J 
            D = diagm(vec(exp.(1/λ * (kron(I(p), E')*v-b-C*x) - ones(p*m, 1)))) # dim: 18x18
            J = [I(p*m)+1/λ*D*C 1/λ*D*kron(I(p), E'); kron(-I(p), E) zeros(p*n, p*n)] # dim: 27x27

            # compute ∇ψ_x
            ∇ψ_x = x - x̂ # 18x1

            # compute ∇̂ψ_b, ∇̂ψ_C
            # ∇̂ψ_b = - 1/λ * [D' zeros(p*m, p*n)] * pinv(J)' * [∇ψ_x; zeros(p*n)]
            # ∇̂ψ_C = - 1/λ * [D' zeros(p*m, p*n)] * pinv(J)' * [∇ψ_x; zeros(p*n)] * x' # 18x1
            z = lsmr(J', [∇ψ_x; zeros(p*n)])
            ∇̂ψ_C = - 1/λ * [D' zeros(p*m, p*n)] * z * x' # 18x1
            push!(∇̂ψ_C_norm, norm(∇̂ψ_C))
            push!(D_norm, norm(D))
            push!(J_norm, norm(J))
            # push!(pinv_J_norm, norm(pinv(J)))

            # update b_plus, C_plus
            # b_plus = proj_B(b - α * ∇̂ψ_b, 0.1)
            C_plus = proj_D(C - α * ∇̂ψ_C, ρ, pa())
        # end
    end

    (;x = x,
      b = b,
      C = C,
      ψ_vals = ψ_vals,
      violation_metrics = violation_metrics,
      lambda_vals = lambda_vals,
      v = v,
      ∇̂ψ_C_norm = ∇̂ψ_C_norm,
      D_norm = D_norm,
      J_norm = J_norm,
      pinv_J_norm = pinv_J_norm,
      F_norm = F_norm)

end

""" RUNNING """
# generate and plot routing games
# function fully_connected_graph()
#     game_name = "fully_connected"

#     edge_dict = Dict(1 => (1,2), 2 => (2,1), 3 => (2,3), 4 => (3,2), 5 => (1,3), 6 => (3,1))
#     g = SimpleDiGraph(3)
#     for i in 1:6; add_edge!(g, edge_dict[i][1], edge_dict[i][2]); end

#     p = 3
#     E = -Matrix(incidence_matrix(g))
#     s = s = vec([1 0 -1; -1 1 0; 0 -1 1])

#     plot_unlabeled(game_name, g)
#     (;game_name=game_name, g=g, p=p, E=E, s=s)
# end

# function simple_graph()
#     game_name = "simple_graph"

#     g = SimpleDiGraph(4)
#     add_edge!(g, 1, 2)
#     add_edge!(g, 1, 3)
#     add_edge!(g, 2, 4)
#     add_edge!(g, 3, 4)
#     add_edge!(g, 2, 3)
#     add_edge!(g, 3, 2)

#     p = 3
#     E = -Matrix(incidence_matrix(g))
#     s = s = vec([1 0 -1; -1 1 0; 0 -1 1])

#     plot_unlabeled(game_name, g)
#     (;game_name=game_name, g=g, p=p, E=E, s=s)
# end

# function grid_graph3_4_players()
#     game_name = "grid_graph3_4_players"

#     row_num = 3
#     col_num = 3
#     g = SimpleDiGraph(row_num * col_num)
#     for i in 1:row_num
#         for j in 1:col_num-1
#         add_edge!(g, row_num*(i-1)+j, row_num*(i-1)+j+1)
#         add_edge!(g, row_num*(i-1)+j+1, row_num*(i-1)+j)
#         end
#     end
#     for j in 1:col_num
#         for i in 1:row_num-1
#             add_edge!(g, row_num*(i-1)+j, row_num*(i-1)+j+row_num)
#             add_edge!(g, row_num*(i-1)+j+row_num, row_num*(i-1)+j)
#         end
#     end

#     p = 4
#     E = -Matrix(incidence_matrix(g))

#     s_p1 = zeros(9); s_p1[4] = 1; s_p1[6] = -1;
#     s_p2 = -s_p1
#     s_p3 = zeros(9); s_p3[2] = 1; s_p3[8] = -1;
#     s_p4 = -s_p3
#     s = [s_p1; s_p2; s_p3; s_p4]

#     p1 = zeros(24); p1[8] = p1[1] = p1[4] = p1[7] = 1;
#     p2 = zeros(24); p2[17] = p2[24] = p2[21] = p2[18] = 1;
#     p3 = zeros(24); p3[4] = p3[7] = p3[17] = p3[24] = 1;
#     p4 = zeros(24); p4[21] = p4[18] = p4[8] = p4[1] = 1;
#     x̂ = vec([p1; p2; p3; p4]) # dim: 24*4 x 1

#     # plot_unlabeled(game_name, g, row_num)
#     (;game_name=game_name, g=g, p=p, E=E, s=s, x̂=x̂)
# end

# function grid_graph5_2_players()
#     game_name = "grid_graph5_2_players"

#     row_num = 5
#     col_num = 5
#     g = SimpleDiGraph(row_num * col_num)
#     for i in 1:row_num
#         for j in 1:col_num-1
#         add_edge!(g, 5*(i-1)+j, 5*(i-1)+j+1)
#         add_edge!(g, 5*(i-1)+j+1, 5*(i-1)+j)
#         end
#     end
#     for j in 1:col_num
#         for i in 1:row_num-1
#             add_edge!(g, 5*(i-1)+j, 5*(i-1)+j+5)
#             add_edge!(g, 5*(i-1)+j+5, 5*(i-1)+j)
#         end
#     end

#     p = 2
#     E = -Matrix(incidence_matrix(g))

#     s_p1 = zeros(25); s_p1[11] = 1; s_p1[15] = -1;
#     s_p2 = -s_p1
#     # s_p3 = zeros(25); s_p3[3] = 1; s_p3[23] = -1;
#     # s_p4 = -s_p3
#     s = [s_p1; s_p2]
#     # s = [s_p1; s_p2; s_p3; s_p4]

#     p1 = zeros(80); p1[32] = p1[15] = p1[19] = p1[23] = p1[27] = p1[31] = 1;
#     p2 = zeros(80); p2[49] = p2[66] = p2[62] = p2[58] = p2[54] = p2[50] = 1;
#     # p3 = zeros(80); p3[7] = p3[11] = p3[28] = p3[46] = p3[64] = p3[77] = 1;
#     # p4 = zeros(80); p4[74] = p4[70] = p4[53] = p4[35] = p4[17] = p4[4] = 1;
#     x̂ = vec([p1; p2]) # dim: 24*2 x 1
#     # x̂ = vec([p1; p2; p3; p4]) # dim: 24*4 x 1

#     # plot_unlabeled(game_name, g)
#     (;game_name=game_name, g=g, p=p, E=E, s=s, x̂=x̂)
# end

# function grid_graph5_4_players()
#     game_name = "grid_graph5_4_players"

#     row_num = 5
#     col_num = 5
#     g = SimpleDiGraph(row_num * col_num)
#     for i in 1:row_num
#         for j in 1:col_num-1
#         add_edge!(g, 5*(i-1)+j, 5*(i-1)+j+1)
#         add_edge!(g, 5*(i-1)+j+1, 5*(i-1)+j)
#         end
#     end
#     for j in 1:col_num
#         for i in 1:row_num-1
#             add_edge!(g, 5*(i-1)+j, 5*(i-1)+j+5)
#             add_edge!(g, 5*(i-1)+j+5, 5*(i-1)+j)
#         end
#     end

#     p = 4
#     E = -Matrix(incidence_matrix(g))

#     s_p1 = zeros(25); s_p1[11] = 1; s_p1[15] = -1;
#     s_p2 = -s_p1
#     s_p3 = zeros(25); s_p3[3] = 1; s_p3[23] = -1;
#     s_p4 = -s_p3
#     s = [s_p1; s_p2; s_p3; s_p4]

#     p1 = zeros(80); p1[32] = p1[15] = p1[19] = p1[23] = p1[27] = p1[31] = 1;
#     p2 = zeros(80); p2[49] = p2[66] = p2[62] = p2[58] = p2[54] = p2[50] = 1;
#     p3 = zeros(80); p3[7] = p3[11] = p3[28] = p3[46] = p3[64] = p3[77] = 1;
#     p4 = zeros(80); p4[74] = p4[70] = p4[53] = p4[35] = p4[17] = p4[4] = 1;
#     x̂ = vec([p1; p2; p3; p4]) # dim: 24*4 x 1

#     # plot_unlabeled(game_name, g)
#     (;game_name=game_name, g=g, p=p, E=E, s=s, x̂=x̂)
# end

# # instantiate a routing game (p, E, s) with desired Nash sol x̂
# game_name, g, p, E, s, x̂ = grid_graph3_4_players()

# assign parameters
# λ = args["lambda"]
# α = args["alpha"]
# ϵ = args["epsilon"]
# ρ = args["rho"]
# max_iter = 10

# # calling method
# println("----------- $(game_name)_λ=($λ)_α=($α)_ϵ=($ϵ)_ρ=($ρ) -----------")
# x, x_init, b, C, ψ_vals = approx_proj_grad(p, E, s, x̂, λ, α, ϵ, ρ, max_iter)

# """ SAVING RESULT """
# # create an individual folder under results/
# dir = "results/$(game_name)/$(game_name)_λ=($λ)_α=($α)_ϵ=($ϵ)_ρ=($ρ)"; mkpath(dir) # mkdir if not exists
# println("saving plot and output to '$dir'")

# # plot and save
# plot(ψ_vals, xlabel="iter", ylabel="ψ(x)"); savefig("$dir/ψ_vals_plots_λ=($λ)_α=($α)_ϵ=($ϵ)_ρ=($ρ).png")
# # @show ψ_vals[end]

# # save outputs
# @save "$dir/output.jld2" x x_init b C ψ_vals
# open("$dir/output.jl", "w") do output_file
#     write(output_file, "# ------------ λ=($λ)_α=($α)_ϵ=($ϵ)_ρ=($ρ) ------------ \n \n")
    
#     write(output_file, "x = ") # writes x =
#     show(output_file, x) # writes the content of x
#     write(output_file, "; \n \n")
    
#     write(output_file, "x_init = ")
#     show(output_file, x_init) 
#     write(output_file, "; \n \n")

#     write(output_file, "b = ")
#     show(output_file, b) 
#     write(output_file, "; \n \n")

#     write(output_file, "C = ")
#     show(output_file, C) 
#     write(output_file, "; \n \n")

#     write(output_file, "ψ_val[end] = ")
#     show(output_file, ψ_vals[end]) 
#     write(output_file, "; \n \n")
# end

# println("--------------------------------\n\n")