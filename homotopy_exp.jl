include("approx_proj_gradient.jl")
include("routing_games.jl")
using ArgParse
using JLD2
using Plots, GraphPlot, ColorSchemes

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
        "--plot"
            help = "plot exp result or not"
            arg_type = Bool
            default = true
        "--forward_solver"
            help = "nlsolve or ipopt"
            arg_type = String
            default = "nlsolve"
    end

    return parse_args(s)
end
args = parse_commandline()


function homotopy_exp_parameter_choice(pa, α, ϵ, ρ, max_iter, λ_list)

    # output placeholder
    x = zeros(pa.p*pa.m)
    v = zeros(pa.p*(pa.n-1))
    b = 0.1*ones(pa.p*pa.m)
    C = zeros(pa.p*pa.m, pa.p*pa.m)
    lambda_vals_list = Vector{Float64}[]
    ψ_vals_exact_list = Float64[]
    ψ_x_exact = Inf

    # run experiment
    println("----------- $(pa.game_name)_ρ=($ρ)_λ_list=($λ_list)_α=($α)_ϵ=($ϵ) -----------\n")
    for λ in λ_list
        if ψ_x_exact <= ϵ
            println("\t----- approx_proj_grad converged -----")
            break
        else
            println("\t----- approx_proj_grad with {λ=$λ} -----\n")
            x, v, b, C, lambda_vals, ψ_vals_exact = approx_proj_grad(pa, λ, α, ϵ, ρ, max_iter, x, v, b, C)
            push!(lambda_vals_list, lambda_vals)
            append!(ψ_vals_exact_list, ψ_vals_exact)
            ψ_x_exact = ψ_vals_exact_list[end]
        end
    end
    println("--------------------------------\n")

     # create an individual folder under homotopy_results/
    dir = "homotopy_results/$(pa.game_name)"
    mkpath(dir) # mkdir if not exists

    # save data to folder
    println("saving result...")
    @save "$dir/$(pa.game_name)_λ_list=($λ_list)_homotopy.jld2" ρ λ_list α ϵ max_iter x v b C lambda_vals_list ψ_vals_exact_list
    println("saved result to '$dir/$(pa.game_name)_λ_list=($λ_list)_homotopy.jld2'")
    println("--------------------------------\n")

    # @load "$dir/$(pa.game_name)_λ_list=($λ_list)_homotopy.jld2" ρ λ_list α ϵ max_iter x v b C lambda_vals_list ψ_vals_exact_list

    # PLOTTING RESULT OR NOT
    if args["plot"] == true
        # plot single axis (log-scaled)
        println("plotting homotopy...")
        plot(title="ρ=($ρ), α=($α), max_iter=($max_iter) \n λ_list=($λ_list)")
        lambda_vals_list_flattened = collect(Iterators.flatten(lambda_vals_list))
        if length(lambda_vals_list) <= 1
            lambda_colors = [:green]
        else
            lambda_colors = palette([:green, :blue], length(lambda_vals_list))
        end
        for i in eachindex(lambda_vals_list)
            lambda_vals = lambda_vals_list[i]
            plot!(findall(x->x==lambda_vals[1], lambda_vals_list_flattened), lambda_vals, yscale=:log10, color = lambda_colors[i], marker=:circle, label="λ=$(lambda_vals[1])", leg=:bottomleft)
        end 
        plot!(ψ_vals_exact_list, yscale=:log10, c=:red, linestyle=:dash, label="ψ(x)_exact", leg=:bottomleft)
        hline!([ϵ], label="ϵ=$ϵ", c=:grey)
        savefig("$dir/$(pa.game_name)_homotopy_log.png")
        println("saved to '$dir/$(pa.game_name)_λ_list=($λ_list)_homotopy_log.png'")

        # plot single axis (linear-scaled)
        plot(title="ρ=($ρ), α=($α), max_iter=($max_iter) \n λ_list=($λ_list)")
        for i in eachindex(lambda_vals_list)
            lambda_vals = lambda_vals_list[i]
            plot!(findall(x->x==lambda_vals[1], lambda_vals_list_flattened), lambda_vals, color = lambda_colors[i], marker=:circle, label="λ=$(lambda_vals[1])", leg=:bottomleft)
        end        
        plot!(ψ_vals_exact_list, c=:red, linestyle=:dash, label="ψ(x)_exact", leg=:bottomleft)
        hline!([ϵ], label="ϵ=$ϵ", c=:grey)
        savefig("$dir/$(pa.game_name)_homotopy_linear.png")
        println("saved to '$dir/$(pa.game_name)_λ_list=($λ_list)_homotopy_linear.png'")
        println("--------------------------------\n")
    end

    (;x = x,
      v = v,
      b = b,
      C = C,
      lambda_vals_list = lambda_vals_list,
      ψ_vals_exact_list = ψ_vals_exact_list)
end


""" RUNNING EXP: PARAMETER CHOICE """
# assign parameters
α = args["alpha"]
ϵ = args["epsilon"]
ρ = args["rho"]
max_iter = args["max_iter"]
λ_list = [1.0, 0.1, 0.01]

# calling method
x, v, b, C, lambda_vals_list, ψ_vals_exact_list = homotopy_exp_parameter_choice(grid_graph3_4_players_reduced(), α, ϵ, ρ, max_iter, λ_list)