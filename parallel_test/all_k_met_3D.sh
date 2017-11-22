#!/bin/bash

/Users/mzumaya/Code/julia/julia -p 4 expansion_folder_par.jl NLOC_MET_3D_EXT 4096 0.01
wait
/Users/mzumaya/Code/julia/julia -p 4 expansion_folder_par.jl NLOC_MET_3D_EXT 4096 0.02
wait
/Users/mzumaya/Code/julia/julia -p 4 expansion_folder_par.jl NLOC_MET_3D_EXT 4096 0.03
wait
/Users/mzumaya/Code/julia/julia -p 4 expansion_folder_par.jl NLOC_MET_3D_EXT 4096 0.25
wait
/Users/mzumaya/Code/julia/julia -p 4 expansion_folder_par.jl NLOC_MET_3D_EXT 4096 0.5
wait
/Users/mzumaya/Code/julia/julia -p 4 expansion_folder_par.jl NLOC_MET_3D_EXT 4096 3.0
wait
/Users/mzumaya/Code/julia/julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.2
wait
/Users/mzumaya/Code/julia/julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.05
wait
/Users/mzumaya/Code/julia/julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.025
wait
/Users/mzumaya/Code/julia/julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.01
wait
/Users/mzumaya/Code/julia/julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.0075
wait
/Users/mzumaya/Code/julia/julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.0035
wait
/Users/mzumaya/Code/julia/julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.0025
wait
/Users/mzumaya/Code/julia/julia -p 4 expansion_folder_par.jl NLOC_TOP_3D_EXT 4096 0.0015
wait
