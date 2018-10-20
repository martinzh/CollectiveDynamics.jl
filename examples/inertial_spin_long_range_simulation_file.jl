### ============== ### ============ ###
##    Extended Inertial Spin Model   ##
##    Martin Zumaya Hernandez        ##
##    EXAMPLE SIMULATION SCRIPT      ##
### ============== ### ============ ###

### ============ INCLUDE PACKAGES ============ ###

using CollectiveDynamics.InertialSpin

### =============== ### =============== ###
###   DEFINITION OF INITIAL PARAMETERS  ###
### =============== ### =============== ###

# N   = 128
# η   = 1.0
# τ   = 6
# rep = 1
# T   = 8*exp10(-5)
# δ   = 0.35

N    = parse(Int, ARGS[1]) # Number of particles
η    = parse(Float64, ARGS[2]) # Dissipation term
T    = parse(Float64, ARGS[3]) # temperature, noise
n_nl = parse(Float64, ARGS[4]) # non local interactions per particle
τ_i  = parse(Int, ARGS[5]) # number of iterations
τ_f  = parse(Int, ARGS[6]) # number of iterations
rep  = parse(Int, ARGS[7]) # repetition number

χ   = 1.25
J   = 0.8
ρ   = 0.3
v0  = 0.1
n_t = 6

pars = InertialParameters(χ, J, η, v0*sqrt(J/χ), ρ, v0, N, T, n_t, n_nl)

# σ = sqrt((2pars.d) * η * T) # noise std deviation ( square root of variance )
σ = sqrt((2pars.d) * η * T * pars.dt) # noise std deviation ( square root of variance )
L = cbrt(N / pars.ρ) # 3D

### =============== ### =============== ###
### SET UP SYSTEM AND OUTPUT STRUCTURE  ###
### =============== ### =============== ###

flock = InertialExtFlock(N, L, pars.v0)

output_path = set_output_data_structure_lr("EXTENDED_INERTIAL_SPIN_EXT", N, η, T, n_nl)

pos_file  = open(output_path * "/pos_$(rep).dat", "w+")
vel_file  = open(output_path * "/vel_$(rep).dat", "w+")
spin_file = open(output_path * "/spin_$(rep).dat", "w+")

### ============== ### ============== ### ============== ###
### INITIALIZATION FROM FILE
### ============== ### ============== ### ============== ###

# pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))
# vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

data_path = joinpath(homedir(),"art_DATA","EXTENDED_INERTIAL_SPIN", "DATA","data_N_$(N)", "eta_$(η)", "eta_$(η)_T_$(T)_nl_$(n_nl)")

raw_data = reinterpret(Float64, read(joinpath(data_path,"pos_$(rep).dat")))
pos_data = raw_data[(end-3N+1):end]

raw_data = reinterpret(Float64, read(joinpath(data_path,"vel_$(rep).dat")))
vel_data = raw_data[(end-3N+1):end]

raw_data = reinterpret(Float64, read(joinpath(data_path,"spin_$(rep).dat")))
spin_data = raw_data[(end-3N+1):end]

for i in 1:N
    flock.pos[i]  = [pos_data[3(i-1)+1], pos_data[3(i-1)+2], pos_data[3(i-1)+3]]
    flock.vel[i]  = [vel_data[3(i-1)+1], vel_data[3(i-1)+2], vel_data[3(i-1)+3]]
    flock.spin[i] = [vel_data[3(i-1)+1], vel_data[3(i-1)+2], vel_data[3(i-1)+3]]
end

# write initial conditions
println("//////// ", 1)
write(pos_file, vcat(flock.pos...))
write(vel_file, vcat(flock.vel...))
write(spin_file, vcat(flock.spin...))

### ============== ### ============== ###
###          SYSTEM EVOLUTION         ###
### ============== ### ============== ###

times = [convert(Int, exp10(i)) for i in τ_i:τ_f]

for i in 1:(length(times) - 1)

    for t in (times[i]+1):times[i+1]

        evolve_extended_system(flock.pos, flock.vel, flock.v_t, flock.v_nl, flock.spin, flock.Rij, pars, σ)

        if t % times[i] == 0 || t % div(times[i], exp10(1)) == 0
            println("//////// ", t)
            write(pos_file, vcat(flock.pos...))
            write(vel_file, vcat(flock.vel...))
            write(spin_file, vcat(flock.spin...))
        end
    end

end

close(pos_file)
close(vel_file)
close(spin_file)

println("Done all")
