#!/bin/bash

for game in "grid_graph3_2_players_reduced" "grid_graph5_4_players_reduced";
do
    for lam in 1.0 0.1 0.01
    do
        julia homotopy_exp.jl --mat true --plot true --game $game --lambda $lam --max_iter 30
    done
done