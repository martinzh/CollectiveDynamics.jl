### ============== ### ============== ### ============== ###
##    Inertial Spin Model Simulation Script               ##
##    Martin Zumaya Hernandez                             ##
##    12 / 03 / 2017                                      ##
### ============== ### ============== ### ============== ###

using CollectiveDynamics

### =============== ### =============== ###
###   DEFINITION OF INITIAL PARAMETERS  ###
### =============== ### =============== ###

# N   = 128
# η   = 1.0
# τ   = 6
# rep = 1
# T   = 8*exp10(-5)

N    = parse(Int, ARGS[1]) # Number of particles
η    = parse(Float64, ARGS[2]) # Dissipation term
T    = parse(Float64, ARGS[3]) # temperature, noise
τ    = parse(Int, ARGS[5]) # number of iterations
rep  = parse(Int, ARGS[6])

χ   = 1.25
J   = 0.8
ρ   = 0.3
v0  = 0.1
δ   = 0.35
n_t = 6

pars = InertialParameters(χ, J, η, v0*sqrt(J/χ), ρ, v0, N, T, n_t, 0)

r0 = (pars.v0 * pars.dt) / pars.l # local interaction range
p  = pars.κ / (N-1) # non-local link probability

σ = sqrt((2pars.d) * η * T) # noise std deviation ( square root of variance )
L = cbrt(N / pars.ρ) # 3D

### =============== ### =============== ###
### SET UP SYSTEM AND OUTPUT STRUCTURE  ###
### =============== ### =============== ###

pos, vel, v_r, v_n, Nij, sp_Nij = set_up_loc_nonLoc_system_3D!(N, L, pars.v0, p)

output_path = set_output_data_structure_3D(N, κ, ω)

pos_file = open(output_path * "/pos_$(rep).dat", "w+")
vel_file = open(output_path * "/vel_$(rep).dat", "w+")
net_file = open(output_path * "/net_$(rep).dat", "w+")

write(net_file, Nij)
close(net_file)

### ============== ### ============== ###
###          SYSTEM EVOLUTION         ###
### ============== ### ============== ###

full_time_evolution_3D(pos_file, vel_file, T, pos, vel, v_r, v_n, sp_Nij, r0, pars.η, ω)

close(pos_file)
close(vel_file)

println("Done all")
