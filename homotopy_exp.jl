include("approx_proj_gradient.jl")
include("routing_games.jl")
using ArgParse
using JLD2
using Plots, GraphPlot

# parse cmd args
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--rho"
            help = "rho"
            arg_type = Float64
            default = 0.5
        "--lambda"
            help = "entropy reg term"
            arg_type = Float64
            default = 0.01
        "--alpha"
            help = "GD step size"
            arg_type = Float64
            default = 0.005
        "--epsilon"
            help = "convergence threshold"
            arg_type = Float64
            default = 1e-3
        "--max_iter"
            help = "maximum iter allowed"
            arg_type = Int64
            default = 10
    end

    return parse_args(s)
end
args = parse_commandline()


function homotopy_exp_parameter_choice(p, E, E_diag, s_reduced, x̂, λ_list, α, ϵ, ρ, max_iter)
    n, m = size(E)

    # output placeholder
    x = zeros(p*m)
    b = 0.1*ones(p*m)
    C = zeros(p*m, p*m)
    v = zeros(p*(n-1))
    
    # ψ_vals_list = Float64[]
    # violation_metrics_list = Float64[]
    lambda_vals_list = Float64[]
    ψ_vals_exact_list = Float64[]


    for λ in λ_list
        println("----- λ=$λ, α=$α -----")
        # x, b, C, ψ_vals, violation_metrics, lambda_vals, v, ∇̂ψ_C_norm, D_norm, J_norm, pinv_J_norm, F_norm, ψ_vals_exact = approx_proj_grad(p, E, E_diag, s_reduced, x̂, λ, α, ϵ, ρ, max_iter, x, b, C)
        x, v, b, C, lambda_vals, ψ_vals_exact = approx_proj_grad(p, E, E_diag, s_reduced, x̂, λ, α, ϵ, ρ, max_iter, x, b, C)
        # append!(ψ_vals_list, ψ_vals)
        # append!(violation_metrics_list, violation_metrics)
        append!(lambda_vals_list, lambda_vals)
        append!(ψ_vals_exact_list, ψ_vals_exact)
    end
    
    (;x = x,
      v = v,
      b = b,
      C = C,
    #   ψ_vals_list = ψ_vals_list,
    #   violation_metrics_list = violation_metrics_list,
      lambda_vals_list = lambda_vals_list,
      ψ_vals_exact_list = ψ_vals_exact_list)
end


""" RUNNING EXP: PARAMETER CHOICE """
# instantiate a routing game (p, E, s) with desired Nash sol x̂
# game_name, g, p, E, E_diag, s_reduced, x̂ = grid_graph3_4_players_reduced()
# game_name, g, p, E, E_diag, s_reduced, x̂ = grid_graph5_2_players_reduced()
game_name, g, p, E, E_diag, s_reduced, x̂ = grid_graph5_4_players_reduced()

# assign parameters
λ = args["lambda"]
α = args["alpha"]
ϵ = args["epsilon"]
ρ = args["rho"]
max_iter = args["max_iter"]
# @show λ_list = [1.0/(2^i) for i in 1:10]
# @show λ_list = [0.002, 0.0016, 0.0015]
# @show λ_list = [0.1, 0.05, 0.025, 0.0125, 0.005, 0.003, 0.002]
@show λ_list = [0.01]

# # calling method
println("----------- $(game_name)_ρ=($ρ)_λ_list=($λ_list)_α=($α)_ϵ=($ϵ) -----------")
x, v, b, C, lambda_vals_list, ψ_vals_exact_list = homotopy_exp_parameter_choice(p, E, E_diag, s_reduced, x̂, λ_list, α, ϵ, ρ, max_iter)


""" SAVING RESULT """
# create an individual folder under homotopy_results/
dir = "homotopy_results/$(game_name)"
mkpath(dir) # mkdir if not exists

# save data to folder
println("saving result to '$dir/$(game_name)_λ_list=($λ_list)_homotopy.jld2'")
@save "$dir/$(game_name)_λ_list=($λ_list)_homotopy.jld2" ρ λ_list α ϵ max_iter x v b C lambda_vals_list ψ_vals_exact_list
println("--------------------------------")


# @load "$dir/$(game_name)_λ_list=($λ_list)_homotopy.jld2" ρ λ_list α ϵ max_iter x v b C lambda_vals_list ψ_vals_exact_list

# plot single axis (log-scaled)
println("plotting homotopy...")
plot(lambda_vals_list, yscale=:log10, c=:green, marker=:circle, label="λ", leg=:bottomleft, xlabel="iter", ylabel="log10 scaled")
# plot!(ψ_vals_list, yscale=:log10, c=:blue, linestyle=:dash, label="ψ(x)", leg=:bottomleft)
# plot!(violation_metrics_list, yscale=:log10, c=:red, linestyle=:dash, label="violation metric", leg=:bottomleft)
plot!(ψ_vals_exact_list, yscale=:log10, c=:purple, linestyle=:dash, label="ψ(x)_exact", leg=:bottomleft)
hline!([ϵ], label="ϵ=$ϵ", c=:grey)
plot!(title="ρ=($ρ), α=($α) \n λ_list=($λ_list)")
savefig("$dir/$(game_name)_λ_list=($λ_list)_homotopy_log.png")
println("saved to '$dir/$(game_name)_λ_list=($λ_list)_homotopy_log.png'")


# plot single axis (linear-scaled)
plot(lambda_vals_list, c=:green, marker=:circle, label="λ", leg=:bottomleft, xlabel="iter", ylabel="linear scaled")
plot!(ψ_vals_exact_list, c=:purple, linestyle=:dash, label="ψ(x)_exact", leg=:bottomleft)
hline!([ϵ], label="ϵ=$ϵ", c=:grey)
plot!(title="ρ=($ρ), α=($α) \n λ_list=($λ_list)")
savefig("$dir/$(game_name)_λ_list=($λ_list)_homotopy_linear.png")
println("saved to '$dir/$(game_name)_λ_list=($λ_list)_homotopy_linear.png'")