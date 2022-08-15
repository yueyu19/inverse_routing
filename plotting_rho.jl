using JLD2
using Plots

rho_list = Float64[]
psi_end_list = Float64[]

for ρ in [1.0, 2.0, 5.0, 10.0]
    @load "results/grid_graph3_λ=(0.01)_α=(0.01)_ϵ=(0.01)_ρ=($ρ)/output.jld2" ψ_vals
    # println("ρ=$ρ", "   ", "ψ_vals[end]=$(ψ_vals[end])")
    push!(rho_list, ρ)
    push!(psi_end_list, ψ_vals[end])
end

println("saving rho paramter search plot in results/")
plot(rho_list, psi_end_list, xlabel="ρ", ylabel="ψ(x) convergence value")
savefig("results/grid_graph3_rho_parameter_search.png")
