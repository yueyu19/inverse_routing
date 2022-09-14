using Graphs, GraphRecipes
using Plots, GraphPlot
using BlockDiagonals

# routing games
# function grid_graph3_2_players()
#     game_name = "grid_graph3_2_players"

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

#     p = 2
#     E = -Matrix(incidence_matrix(g))
#     n, m = size(E)

#     s_p1 = zeros(9); s_p1[4] = 1; s_p1[6] = -1;
#     s_p2 = -s_p1
#     s = [s_p1; s_p2]

#     p1 = zeros(24); p1[8] = p1[1] = p1[4] = p1[7] = 1;
#     p2 = zeros(24); p2[17] = p2[24] = p2[21] = p2[18] = 1;
#     x̂ = vec([p1; p2]) # dim: 24*2

#     # plot_unlabeled(game_name, g, row_num)
#     (;game_name=game_name, g=g, p=p, E=E, n=n, m=m, s=s, x̂=x̂)
# end

function grid_graph3_2_players_reduced()
    game_name = "grid_graph3_2_players_reduced"

    row_num = 3
    col_num = 3
    g = SimpleDiGraph(row_num * col_num)
    for i in 1:row_num
        for j in 1:col_num-1
        add_edge!(g, row_num*(i-1)+j, row_num*(i-1)+j+1)
        add_edge!(g, row_num*(i-1)+j+1, row_num*(i-1)+j)
        end
    end
    for j in 1:col_num
        for i in 1:row_num-1
            add_edge!(g, row_num*(i-1)+j, row_num*(i-1)+j+row_num)
            add_edge!(g, row_num*(i-1)+j+row_num, row_num*(i-1)+j)
        end
    end

    p = 2

    E = -Matrix(incidence_matrix(g)); n, m = size(E)
    E1 = E[1:end .∉ [[6]], 1:end]
    E2 = E[1:end .∉ [[4]], 1:end]
    E_diag = BlockDiagonal([E1, E2])

    n, m = size(E)

    s_p1 = zeros(9); s_p1[4] = 1; s_p1[6] = -1;
    s_p2 = -s_p1
    s = [s_p1; s_p2]
    s_reduced = filter!(i->i!=-1, s)

    p1 = zeros(24); p1[8] = p1[1] = p1[4] = p1[7] = 1;
    p2 = zeros(24); p2[17] = p2[24] = p2[21] = p2[18] = 1;
    x̂ = vec([p1; p2]) # dim: 24*2

    # plot_unlabeled(game_name, g, row_num)
    (;game_name=game_name, g=g, p=p, E=E, n=n, m=m, E_diag=E_diag, s_reduced=s_reduced, x̂=x̂)
end

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
#     n, m = size(E)

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
#     (;game_name=game_name, g=g, p=p, E=E, n=n, m=m, s=s, x̂=x̂)
# end

# function grid_graph3_4_players_reduced()
#     game_name = "grid_graph3_4_players_reduced"

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

#     E = -Matrix(incidence_matrix(g)); n, m = size(E)
#     E1 = E[1:end .∉ [[6]], 1:end]
#     E2 = E[1:end .∉ [[4]], 1:end]
#     E3 = E[1:end .∉ [[8]], 1:end]
#     E4 = E[1:end .∉ [[2]], 1:end]
#     E_diag = BlockDiagonal([E1, E2, E3, E4])
    
#     s_p1 = zeros(9); s_p1[4] = 1; s_p1[6] = -1;
#     s_p2 = -s_p1
#     s_p3 = zeros(9); s_p3[2] = 1; s_p3[8] = -1;
#     s_p4 = -s_p3
#     s = [s_p1; s_p2; s_p3; s_p4]
#     s_reduced = filter!(i->i!=-1, s)

#     p1 = zeros(24); p1[8] = p1[1] = p1[4] = p1[7] = 1;
#     p2 = zeros(24); p2[17] = p2[24] = p2[21] = p2[18] = 1;
#     p3 = zeros(24); p3[4] = p3[7] = p3[17] = p3[24] = 1;
#     p4 = zeros(24); p4[21] = p4[18] = p4[8] = p4[1] = 1;
#     x̂ = vec([p1; p2; p3; p4]) # dim: 24*4 x 1

#     # plot_unlabeled(game_name, g, row_num)
#     (;game_name=game_name, g=g, p=p, E=E, n=n, m=m, E_diag=E_diag, s_reduced=s_reduced, x̂=x̂)
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
#     E = -Matrix(incidence_matrix(g)); n, m = size(E)

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
#     (;game_name=game_name, g=g, p=p, E=E, n=n, m=m, s=s, x̂=x̂)
# end

# function grid_graph5_2_players_reduced()
#     game_name = "grid_graph5_2_players_reduced"

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
#     E = -Matrix(incidence_matrix(g)); n, m = size(E)
#     E1 = E[1:end .∉ [[15]], 1:end]
#     E2 = E[1:end .∉ [[11]], 1:end]
#     E_diag = BlockDiagonal([E1, E2])

#     s_p1 = zeros(25); s_p1[11] = 1; s_p1[15] = -1;
#     s_p2 = -s_p1
#     # s_p3 = zeros(25); s_p3[3] = 1; s_p3[23] = -1;
#     # s_p4 = -s_p3
#     s = [s_p1; s_p2]
#     # s = [s_p1; s_p2; s_p3; s_p4]
#     s_reduced = filter!(i->i!=-1, s)
    
#     p1 = zeros(80); p1[32] = p1[15] = p1[19] = p1[23] = p1[27] = p1[31] = 1;
#     p2 = zeros(80); p2[49] = p2[66] = p2[62] = p2[58] = p2[54] = p2[50] = 1;
#     # p3 = zeros(80); p3[7] = p3[11] = p3[28] = p3[46] = p3[64] = p3[77] = 1;
#     # p4 = zeros(80); p4[74] = p4[70] = p4[53] = p4[35] = p4[17] = p4[4] = 1;
#     x̂ = vec([p1; p2]) # dim: 24*2 x 1
#     # x̂ = vec([p1; p2; p3; p4]) # dim: 24*4 x 1

#     # plot_unlabeled(game_name, g)
#     (;game_name=game_name, g=g, p=p, E=E, n=n, m=m, E_diag=E_diag, s_reduced=s_reduced, x̂=x̂)
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
#     E = -Matrix(incidence_matrix(g)); n, m = size(E)

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
#     (;game_name=game_name, g=g, p=p, E=E, n=n, m=m, s=s, x̂=x̂)
# end

function grid_graph5_4_players_reduced()
    game_name = "grid_graph5_4_players_reduced"

    row_num = 5
    col_num = 5
    g = SimpleDiGraph(row_num * col_num)
    for i in 1:row_num
        for j in 1:col_num-1
        add_edge!(g, 5*(i-1)+j, 5*(i-1)+j+1)
        add_edge!(g, 5*(i-1)+j+1, 5*(i-1)+j)
        end
    end
    for j in 1:col_num
        for i in 1:row_num-1
            add_edge!(g, 5*(i-1)+j, 5*(i-1)+j+5)
            add_edge!(g, 5*(i-1)+j+5, 5*(i-1)+j)
        end
    end

    p = 4
    E = -Matrix(incidence_matrix(g)); n, m = size(E)
    E1 = E[1:end .∉ [[15]], 1:end]
    E2 = E[1:end .∉ [[11]], 1:end]
    E3 = E[1:end .∉ [[18]], 1:end]
    E4 = E[1:end .∉ [[8]], 1:end]
    E_diag = BlockDiagonal([E1, E2, E3, E4])

    s_p1 = zeros(25); s_p1[11] = 1; s_p1[15] = -1;
    s_p2 = -s_p1
    s_p3 = zeros(25); s_p3[8] = 1; s_p3[18] = -1;
    s_p4 = -s_p3
    s = [s_p1; s_p2; s_p3; s_p4]
    s_reduced = filter!(i->i!=-1, s)

    p1 = zeros(80); p1[32] = p1[14] = p1[1] = p1[4] = p1[7] = p1[10] = p1[13] = p1[31] = 1;
    p2 = zeros(80); p2[49] = p2[67] = p2[80] = p2[77] = p2[74] = p2[71] = p2[68] = p2[50] = 1;
    p3 = zeros(80); p3[23] = p3[28] = p3[46] = p3[62] = 1;
    p4 = zeros(80); p4[58] = p4[53] = p4[35] = p4[19] = 1;
    x̂ = vec([p1; p2; p3; p4]) # dim: 24*4 x 1

    # plot_unlabeled(game_name, g, row_num)
    (;game_name=game_name, g=g, p=p, E=E, n=n, m=m, E_diag=E_diag, s_reduced=s_reduced, x̂=x̂)
end

# # helper function: plot unlabeled directed graph - GRID_GRAPH
function plot_unlabeled(game_name, g, row_num)
    dir = "results/$(game_name)"
    mkpath(dir)

    # generating edgelabel_dict
    edgelabel_dict = Dict{Tuple{Int64, Int64}, String}()
    for i in 1:ne(g)
        col = -Matrix(incidence_matrix(g))[:,i]
        edgelabel_dict[(findall(x->x==1, col)[1], findall(x->x==-1, col)[1])] =string(i)
    end
    # plotting directed graph with node and link labelings
    δ=0.01
    locs_x = collect(Iterators.flatten([i%2==1 ? [1:row_num...] : [1:row_num...].+δ for i in 1:row_num]))
    locs_y = collect(Iterators.flatten([repeat([i],row_num) for i in row_num:-1:1]))
    graphplot(g,
              names=1:nv(g), x=2*locs_x, y=1.5*locs_y, nodeshape=:circle, nodesize=0.8,
              edgelabel=edgelabel_dict, curvature_scalar=0.2)
    savefig("$dir/($game_name$row_num)_unlabeled.png")
end

# # helper function: plot labeled directed graph
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
