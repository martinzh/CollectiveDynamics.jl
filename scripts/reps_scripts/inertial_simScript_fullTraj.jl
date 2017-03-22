### ============== ### ============== ### ============== ###
##    Inertial Spin Model Simulation Script               ##
##    Martin Zumaya Hernandez                             ##
##    12 / 03 / 2017                                      ##
### ============== ### ============== ### ============== ###

using CollectiveDynamics, CollectiveDynamics.DataAnalysis

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
τ    = parse(Int, ARGS[4]) # number of iterations
rep  = parse(Int, ARGS[5]) # repetition number

χ   = 1.25
J   = 0.8
ρ   = 0.3
v0  = 0.1
n_t = 6

# δ   = 0.35

pars = InertialParameters(χ, J, η, v0*sqrt(J/χ), ρ, v0, N, T, n_t, 0)

σ = sqrt((2pars.d) * η * T) # noise std deviation ( square root of variance )
L = cbrt(N / pars.ρ) # 3D

### =============== ### =============== ###
### SET UP SYSTEM AND OUTPUT STRUCTURE  ###
### =============== ### =============== ###

flock = InertialFlock(N, L, pars.v0)

make_dir_from_path("$(homedir())/DATA_TEST")
make_dir_from_path("$(homedir())/DATA_TEST/non_loc_N_$(N)_eta_$(η)_T_$(T))")

output_path = "$(homedir())/DATA_TEST/non_loc_N_$(N)_eta_$(η)_T_$(T))"

pos_file  = open(output_path * "/pos_$(rep).dat", "w+")
vel_file  = open(output_path * "/vel_$(rep).dat", "w+")
spin_file = open(output_path * "/spin_$(rep).dat", "w+")

### ============== ### ============== ###
###          SYSTEM EVOLUTION         ###
### ============== ### ============== ###

# full_time_evolution_inertial_system(pos_file, vel_file, spin_file, τ, flock, pars, σ)

for t in 1:convert(Int, exp10(τ))

    CollectiveDynamics.evolve_inertial_system(flock.pos, flock.vel, flock.v_t, flock.spin, flock.Rij, pars, σ)

    println("//////// ", t)
    write(pos_file, vcat(flock.pos...))
    write(vel_file, vcat(flock.vel...))
    write(spin_file, vcat(flock.spin...))
end

close(pos_file)
close(vel_file)
close(spin_file)

println("Done all")
