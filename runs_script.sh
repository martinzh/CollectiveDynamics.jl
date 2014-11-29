#!/bin/bash
nohup julia Defs.jl 0.05 &
nohup julia Defs.jl 0.1 &
nohup julia Defs.jl 0.15 &
julia Defs.jl 0.2
nohup julia Defs.jl 0.25 &
nohup julia Defs.jl 0.3 &
nohup julia Defs.jl 0.35 &
julia Defs.jl 0.4
nohup julia Defs.jl 0.45 &
nohup julia Defs.jl 0.5 &
nohup julia Defs.jl 0.55 &
julia Defs.jl 0.6
