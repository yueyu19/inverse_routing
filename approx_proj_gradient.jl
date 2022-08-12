using LinearAlgebra
include("entropy_routing_solver.jl")
using Plots
using Graphs, GraphRecipes
using Colors
using ArgParse

# parse cmd args
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--rho"
            help = "rho"
            arg_type = Float64
            default = 1.0
    end

    return parse_args(s)
end
args = parse_commandline()


# helper function: projection onto D
function proj_D(C, ρ)
    s = eigvals(0.5 * (C + C'))
    U = eigvecs(0.5 * (C + C'))
    A = 1/2 * (C - C') + U * diagm(vec(max.(s,zeros(length(s))))) * U'
    return ρ / max.(ρ, norm(A)) * A
end

"""
Approximate projected gradient method
Inputs:
- ψ, measuring function, e.g. ψ(x)=1/2*||x-x̂||^2
- λ, entropy weight
- α, step size
- ϵ, stopping tolerance
- ρ

Returns:
- x: mixed eq strategies for each player on space of links
- b
- C
"""
function approx_proj_grad(p, E, s, λ, α, ϵ, ρ, max_iter)
    n, m = size(E)
    x̂ = [0, 0, 0, 1, 1, 0, 
        0, 1, 0, 0, 1, 0, 
        0, 1, 0, 1, 0, 0] # 18x1
    
    # initiate b, b_plus, C, C_plus
    b = zeros(p*m)
    b_plus = 0.1*ones(p*m)
    C = 0.1*I(p*m) 
    C_plus = zeros(p*m, p*m)

    ψ_vals = Float64[]
    x_init = Vector{Float64}
    x = Vector{Float64}

    println("Starting approx proj grad...")

    for i in 1:max_iter
        if max(norm(b-b_plus), norm(C-C_plus)) <= ϵ
            break
        else
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
            if i == 1; x_init = x; end # store initial x
            push!(ψ_vals, 0.5 * norm(x-x̂)^2)

            # compute D, J 
            D = diagm(vec(exp.(1/λ * (kron(I(p), E')*v-b-C*x) - ones(p*m, 1)))) # dim: 18x18
            J = [I(p*m)+1/λ*D*C 1/λ*D*kron(I(p), E'); kron(-I(p), E) zeros(p*n, p*n)] # dim: 27x27

            # compute ∇ψ_x
            ∇ψ_x = x - x̂ # 18x1

            # compute ∇̂ψ_b, ∇̂ψ_C
            # ∇̂ψ_b = - 1/λ * [D' zeros(p*m, p*n)] * pinv(J)' * [∇ψ_x; zeros(p*n)]
            ∇̂ψ_C = - 1/λ * [D' zeros(p*m, p*n)] * pinv(J)' * [∇ψ_x; zeros(p*n)] * x' # 18x1

            # update b_plus, C_plus
            # b_plus = proj_B(b - α * ∇̂ψ_b, 0.1)
            C_plus = proj_D(C - α * ∇̂ψ_C, ρ)
        end
    end
    println("Finished.")

    (;x = x,
      x_init = x_init,
      b = b,
      C = C,
      ψ_vals = ψ_vals)

end

# create routing game instance directly 
# p = 3
# E = [1 -1  0  0  1 -1; -1  1  1 -1  0  0; 0  0 -1  1 -1 1]
# @show E
# s = vec([1 0 -1; -1 1 0; 0 -1 1])

# # create routing network
game_name = "fully_connected"
edge_dict = Dict(1 => (1,2), 2 => (2,1), 3 => (2,3), 4 => (3,2), 5 => (1,3), 6 => (3,1))
g = SimpleDiGraph(3)
for i in 1:6
    add_edge!(g, edge_dict[i][1], edge_dict[i][2])
end

# create routing game instance
p = 3
E = -Matrix(incidence_matrix(g))
s = s = vec([1 0 -1; -1 1 0; 0 -1 1])

# assign parameters
λ = 0.01
α = 0.01
ϵ = 0.01
ρ = args["rho"]
max_iter = 10

# calling method
# x, x_init, b, C, ψ_vals = approx_proj_grad(p, E, s, λ, α, ϵ, ρ, max_iter)
# plot(ψ_vals, xlabel="iter", ylabel="ψ(x)")
# savefig("ψ_vals_plots_λ=($λ)_α=($α)_ϵ=($ϵ)_ρ=($ρ).png")
# @show ψ_vals[end]

# plot unlabeled
function plot_unlabeled(game_name, g)
    # generating edgelabel_dict
    edgelabel_dict = Dict{Tuple{Int64, Int64}, String}()
    for i in 1:ne(g)
        col = -Matrix(incidence_matrix(g))[:,i]
        edgelabel_dict[(findall(x->x==1, col)[1], findall(x->x==-1, col)[1])] =string(i)
    end
    # plotting directed graph with node and link labelings
    graphplot(g, fontsize=11,
              names=1:nv(g), nodeshape=:circle, nodesize=0.2,
              edgelabel=edgelabel_dict)
    savefig("($game_name)_unlabeled.png")
end

# plot labeled
function plot_labeled(game_name, g, x, player)
    # highlight referred player
    node_colors = [colorant"grey" for i in 1:nv(g)]; node_colors[player] = colorant"red"
    # highlight chosen path from x
    edgelabel_dict = Dict((1,2) => "1", (2,1) => "2", (2,3) => "3", (3,2) => "4", (1,3) => "5", (3,1) => "6")
    for i in findall(x->x>0.99, x[(player-1)*ne(g)+1 : player*ne(g)])
        edgelabel_dict[edge_dict[i]] *= "*"
    end
    @show x[(player-1)*ne(g)+1 : player*ne(g)]
    # plot labeled graph
    graphplot(g, fontsize=11,
            names=1:nv(g), nodeshape=:circle, nodesize=0.2, markercolor = node_colors,
            edgelabel=edgelabel_dict)
    savefig("($game_name)_labeled_player_$player.png")
end

plot_unlabeled(game_name, g)
# plot_labeled(game_name*"_x_init", g, x_init, 2)
# plot_labeled(game_name*"_x", g, x, 2)


#todo: 5x5 grid world like directed graph, p=4, s=collision source-sink Vector