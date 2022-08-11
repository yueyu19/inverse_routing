#!/bin/bash

for r in 0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
    echo "rho = $r:"
    julia approx_proj_gradient.jl --rho $r
    echo "--------------------------------"
done