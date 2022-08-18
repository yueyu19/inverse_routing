using LinearAlgebra
include("entropy_routing_solver.jl")
using Graphs, GraphRecipes

function grid_graph3()
    game_name = "grid_graph"

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
    @show nv(g)
    @show ne(g)

    p = 4
    E = -Matrix(incidence_matrix(g))

    p1 = zeros(9); p1[4] = 1; p1[6] = -1;
    p2 = -p1
    p3 = zeros(9); p3[2] = 1; p3[8] = -1;
    p4 = -p3
    s = [p1;p2;p3;p4]

    # plot_unlabeled(game_name, g, row_num)
    (;game_name=game_name, g=g, p=p, E=E, s=s)
end

function grid_graph5()
    game_name = "grid_graph"

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
    @show nv(g)
    @show ne(g)

    E = -Matrix(incidence_matrix(g))

    p = 4
    p1 = zeros(25); p1[11] = 1; p1[15] = -1;
    p2 = -p1
    p3 = zeros(25); p3[3] = 1; p3[23] = -1;
    p4 = -p3
    s = [p1;p2;p3;p4]

    # plot_unlabeled(game_name, g)
    (;game_name=game_name, g=g, p=p, E=E, s=s)
end

game_name, g, p, E, s = grid_graph5()
λ = 0.01

n, m = size(E)
b = 0.1*ones(p*m)
C = zeros(p*m, p*m)

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
@time x, v = solve_entropy_routing(pa(), λ)
# @show x
# @show v

# F = [x - exp.(1/λ * (kron(I(p), E') * v - b - C * x) - ones(p*m,1));
#     s - kron(I(p), E) * x]
# @show size(F)
# @show max(abs.(F)...)