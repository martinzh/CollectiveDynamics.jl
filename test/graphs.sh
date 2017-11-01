#!/bin/bash

F="$1"
k="$2"

nohup ~/source/julia-0.6.0/julia graph_grid_trays.jl $F 512 $k 0.025 &
nohup ~/source/julia-0.6.0/julia graph_grid_trays.jl $F 512 $k 0.05 &
nohup ~/source/julia-0.6.0/julia graph_grid_trays.jl $F 512 $k 0.075 &
nohup ~/source/julia-0.6.0/julia graph_grid_trays.jl $F 512 $k 0.1 &
nohup ~/source/julia-0.6.0/julia graph_grid_trays.jl $F 512 $k 0.25 &
nohup ~/source/julia-0.6.0/julia graph_grid_trays.jl $F 512 $k 0.6 &
nohup ~/source/julia-0.6.0/julia graph_grid_trays.jl $F 512 $k 0.65 &
nohup ~/source/julia-0.6.0/julia graph_grid_trays.jl $F 512 $k 0.7 &
nohup ~/source/julia-0.6.0/julia graph_grid_trays.jl $F 512 $k 0.8 &
nohup ~/source/julia-0.6.0/julia graph_grid_trays.jl $F 512 $k 1.0 &
