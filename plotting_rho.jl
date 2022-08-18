using JLD2
using Plots

game_name = "grid_graph5"

# parse rho_list from cmd line arg
rho_list = let expr = Meta.parse(ARGS[1])
    @assert expr.head == :vect
    Float64.(expr.args)
end

# collect psi_end_list
psi_end_list = Float64[]
for ρ in rho_list
    @load "results/$(game_name)/$(game_name)_λ=(0.01)_α=(0.01)_ϵ=(0.01)_ρ=($ρ)/output.jld2" ψ_vals
    println("ρ=$ρ", "   ", "ψ_vals[end]=$(ψ_vals[end])")
    push!(psi_end_list, ψ_vals[end])
end

# plot and save rho graph in results/
println("saving rho paramter search plot in 'results/$(game_name)_rho_parameter_search.png'")
plot(rho_list, psi_end_list, xlabel="ρ", ylabel="ψ(x) convergence value")
savefig("results/$(game_name)/$(game_name)_rho_parameter_search.png")
