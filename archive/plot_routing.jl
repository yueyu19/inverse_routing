using Graphs
using GraphRecipes, Plots

# function graphToMatrix(g)
#     im = Matrix(incidence_matrix(g))
#     am = Matrix(adjacency_matrix(g))
#     (;im=im, am=am)
# end

# create directed graph
g = SimpleDiGraph(4)
add_edge!(g, 1, 2)
add_edge!(g, 1, 3)
add_edge!(g, 2, 4)
add_edge!(g, 3, 4)
add_edge!(g, 2, 3)
add_edge!(g, 3, 2)

# create incidence matrix from directed graph for computation
im = Matrix(incidence_matrix(g))
E = -im; @show E

# plot directed graph
δ=0.01
locs_x = [1, 2-δ, 2+δ, 3]
locs_y = [2, 1, 3, 2]
edgelabel_dict = Dict((1,2) => "1", (1,3) => "2", (2,4) => "3", (3,4) => "4", (2,3) => "5", (3,2) => "6")
graphplot(g, fontsize=11,
          names=1:nv(g), x=locs_x, y=locs_y, nodeshape=:circle, nodesize=0.4,
          edgelabel=edgelabel_dict, curvature_scalar=0.1)
savefig("routing_plot_coordinated.png")
