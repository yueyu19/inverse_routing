using JLD2, DelimitedFiles
using CSV, DataFrames, MAT
include("routing_solvers/routing_solver_exact_jump.jl")
include("routing_games.jl")

"""lambda plots"""
# 3x3_2p
# @load "results/grid_graph3_2_players_reduced/grid_graph3_2_players_reduced_λ_list=([1.0])_homotopy.jld2" ρ λ_list α ϵ max_iter x v b C lambda_vals_list ψ_vals_exact_list
# ψ_vals_exact_list_lambda_1 = ψ_vals_exact_list
# @load "results/grid_graph3_2_players_reduced/grid_graph3_2_players_reduced_λ_list=([0.1])_homotopy.jld2" ρ λ_list α ϵ max_iter x v b C lambda_vals_list ψ_vals_exact_list
# ψ_vals_exact_list_lambda_0_1 = ψ_vals_exact_list
# @load "results/grid_graph3_2_players_reduced/grid_graph3_2_players_reduced_λ_list=([0.01])_homotopy.jld2" ρ λ_list α ϵ max_iter x v b C lambda_vals_list ψ_vals_exact_list
# ψ_vals_exact_list_lambda_0_0_1 = ψ_vals_exact_list

# df = DataFrame(lambda_1 = ψ_vals_exact_list_lambda_1)
#             #    lambda_0_1 = ψ_vals_exact_list_lambda_0_1,
#             #    lambda_0_0_1 = ψ_vals_exact_list_lambda_0_0_1)
# CSV.write("3x3_2p_lambda_plot.csv", df)

# 5x5_4p
# @load "results/grid_graph5_4_players_reduced/grid_graph5_4_players_reduced_λ_list=([1.0])_homotopy.jld2" ρ λ_list α ϵ max_iter x v b C lambda_vals_list ψ_vals_exact_list
# ψ_vals_exact_list_lambda_1 = ψ_vals_exact_list
# @load "results/grid_graph5_4_players_reduced/grid_graph5_4_players_reduced_λ_list=([0.1])_homotopy.jld2" ρ λ_list α ϵ max_iter x v b C lambda_vals_list ψ_vals_exact_list
# ψ_vals_exact_list_lambda_0_1 = ψ_vals_exact_list
# @load "results/grid_graph5_4_players_reduced/grid_graph5_4_players_reduced_λ_list=([0.01])_homotopy.jld2" ρ λ_list α ϵ max_iter x v b C lambda_vals_list ψ_vals_exact_list
# ψ_vals_exact_list_lambda_0_0_1 = ψ_vals_exact_list

# df1 = DataFrame(lambda_1 = ψ_vals_exact_list_lambda_1)
# df2 = DataFrame(lambda_0_1 = ψ_vals_exact_list_lambda_0_1)
# df3 = DataFrame(lambda_0_0_1 = ψ_vals_exact_list_lambda_0_0_1)
# CSV.write("5x5_4p_lambda_plot.csv", df3)

"""rho plots"""
# # 5x5 4p
# rho_list = [0.01, 0.05, 0.1, 0.3, 0.5]
# # 0.08, 0.1, 0.2, 0.3, 0.5, 0.7] 

# alpha_list = 0.01 * ones(length(rho_list));

# epsilon_list = 0.01 * ones(length(rho_list)); 
# epsilon_list[findall(x->x==0.1, rho_list)[1]] = 0.1; 
# epsilon_list[findall(x->x==0.5, rho_list)[1]] = 0.1; 
# epsilon_list[findall(x->x==0.3, rho_list)[1]] = 0.1; 

# game_name = "grid_graph5_4_players"
# λ = 0.01
# psi_end_list = Float64[]
# for i in 1:length(rho_list)
#     ρ = rho_list[i]
#     α = alpha_list[i]
#     ϵ = epsilon_list[i]
#     @load "rho_exp/results/$(game_name)/$(game_name)_λ=($λ)_α=($α)_ϵ=($ϵ)_ρ=($ρ)/output.jld2" ψ_vals
#     println("ρ=$ρ", "   ", "ψ_vals[end]=$(ψ_vals[end])")
#     push!(psi_end_list, ψ_vals[end])
# end

# df = DataFrame(rho = rho_list,
#                 psi_convergence_val = psi_end_list)
# CSV.write("5x5_4p_rho_plot.csv", df)

"""path cost bar plots"""



"""save exact x, v to file too"""
pa = grid_graph5_4_players_reduced()
λ_list = [1.0]
@load "results/$(pa.game_name)/$(pa.game_name)_λ_list=($λ_list)_homotopy.jld2" ρ λ_list α ϵ max_iter x v b C lambda_vals_list ψ_vals_exact_list
x_exact, v_exact = solve_routing(pa, b, C, zeros(pa.p*pa.m), zeros(pa.p*(pa.n-1)))
@save "results/$(pa.game_name)/$(pa.game_name)_λ_list=($λ_list)_homotopy.jld2" ρ λ_list α ϵ max_iter x v b C lambda_vals_list ψ_vals_exact_list x_exact v_exact

file = matopen("results/$(pa.game_name)/$(pa.game_name)_λ_list=($λ_list)_homotopy.mat", "w")
write(file, "lambda_vals_list", lambda_vals_list)
write(file, "lambda_vals_list_flattened", collect(Iterators.flatten(lambda_vals_list)))
write(file, "psi_vals_exact_list", ψ_vals_exact_list)
write(file, "x", x)
write(file, "v", v)
write(file, "b", b)
write(file, "C", C)
write(file, "x_exact", x_exact)
write(file, "v_exact", v_exact)
close(file)