#!/bin/bash

julia -p 4 expansion_folder_par.jl NLOC_MET_3D_EXT 4096 0.01
wait
julia -p 4 expansion_folder_par.jl NLOC_MET_3D_EXT 4096 0.02
wait
julia -p 4 expansion_folder_par.jl NLOC_MET_3D_EXT 4096 0.03
wait
julia -p 4 expansion_folder_par.jl NLOC_MET_3D_EXT 4096 0.25
wait
julia -p 4 expansion_folder_par.jl NLOC_MET_3D_EXT 4096 0.5
wait
julia -p 4 expansion_folder_par.jl NLOC_MET_3D_EXT 4096 3.0
wait
julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.2
wait
julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.05
wait
julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.025
wait
julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.01
wait
julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.0075
wait
julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.0035
wait
julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.0025
wait
julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.0015
wait
