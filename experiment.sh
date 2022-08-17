#!/bin/bash

# for r in 0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
# for r in 0.005 0.01 0.05 0.1 0.5 1 3 5 10
for r in 0.01 0.05 0.1 0.5 1 5 10
do
    # echo "rho = $r:"
    julia approx_proj_gradient.jl --rho $r
    # echo "--------------------------------\n\n"
done

julia plotting_rho.jl "[0.01, 0.05, 0.1, 0.5, 1, 5, 10]"