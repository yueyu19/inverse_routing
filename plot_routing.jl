using LightGraphs
using GraphRecipes, Plots

# create directed graph
g = SimpleDiGraph(4)
add_edge!(g, 1, 2)
add_edge!(g, 1, 3)
add_edge!(g, 2, 4)
add_edge!(g, 3, 4)
add_edge!(g, 2, 3)
add_edge!(g, 3, 2)

function graphToMatrix(g)
    im = Matrix(incidence_matrix(g))
    am = Matrix(adjacency_matrix(g))
    (;im=im, am=am)
end

# create incidence matrix from directed graph for computation
im = Matrix(incidence_matrix(g))
# E = [1 1 0 0 0 0; -1 0 1 0 1 -1; 0 -1 0 1 -1 1; 0 0 -1 -1 0 0]
E = -im
@show E

# create adjacency matrix from directed graph for visualization
am = Matrix(adjacency_matrix(g))

# plot directed graph
graphplot(am, names=1:4)
savefig("routing_plot.png")
