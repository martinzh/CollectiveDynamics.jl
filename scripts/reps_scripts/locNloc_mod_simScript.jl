### ============== ### ============== ### ============== ###
##    2D Local and NonLocal Model Simulation Script       ##
##    Martin Zumaya Hernandez                             ##
##    12 / 03 / 2017                                      ##
### ============== ### ============== ### ============== ###

using CollectiveDynamics

### =============== ### =============== ###
###   DEFINITION OF INITIAL PARAMETERS  ###
### =============== ### =============== ###

# N = 128
# κ = 2.0
# ω = 0.5
# T = 3
# rep = 1

N   = parse(Int, ARGS[1]) # number of particles
κ   = parse(Float64, ARGS[2]) # average non-local interactions
ω   = parse(Float64, ARGS[3]) # interactions relative weight
T   = parse(Int, ARGS[4]) # integration time steps

dim = parse(Int, ARGS[5])
mod = ARGS[6]

rep = parse(Int, ARGS[7])

ρ = 0.3
η = 0.15

pars = LocNonLocParameters(N, κ, ω, ρ, η) # system's parameters

L  = sqrt(N / pars.ρ) # size of box
r0 = (pars.v0 * pars.dt) / pars.l # local interaction range
p  = pars.κ / (N-1) # non-local link probability

### =============== ### =============== ###
### SET UP SYSTEM AND OUTPUT STRUCTURE  ###
### =============== ### =============== ###

flock = LocNonLocFlock(N, L, pars.v0, p, dim)

if dim == 2 && mod == "mod"
    output_path = set_output_data_structure_2D_MOD(N, κ, ω)
elseif dim == 2 && mod == "org"
    output_path = set_output_data_structure_2D(N, κ, ω)
elseif dim == 3 && mod == "mod"
    output_path = set_output_data_structure_3D_MOD(N, κ, ω)
elseif dim == 3 && mod == "org"
    output_path = set_output_data_structure_3D(N, κ, ω)
end

pos_file = open(output_path * "/pos_$(rep).dat", "w+")
vel_file = open(output_path * "/vel_$(rep).dat", "w+")
net_file = open(output_path * "/net_$(rep).dat", "w+")

write(net_file, flock.Nij)
close(net_file)

### ============== ### ============== ###
###          SYSTEM EVOLUTION         ###
### ============== ### ============== ###


if dim == 2 && mod == "mod"
    full_time_evolution_2D_MOD(pos_file, vel_file, T, flock, r0, pars.η, ω)
elseif dim == 2 && mod == "org"
    full_time_evolution_2D(pos_file, vel_file, T, flock, r0, pars.η, ω)
elseif dim == 3 && mod == "mod"
    full_time_evolution_3D_MOD(pos_file, vel_file, T, flock, r0, pars.η, ω)
elseif dim == 3 && mod == "org"
    full_time_evolution_3D(pos_file, vel_file, T, flock, r0, pars.η, ω)
end

close(pos_file)
close(vel_file)

println("Done all")
