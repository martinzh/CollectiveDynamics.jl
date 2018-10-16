### ============== ### ============ ###
##    Extended Inertial Spin Model   ##
##    (SHARED ARRAYS)                ##
##    Martin Zumaya Hernandez        ##
##    EXAMPLE SIMULATION SCRIPT      ##
### ============== ### ============ ###

### ============ INCLUDE PACKAGES ============ ###

@everywhere using CollectiveDynamics.InertialSpin

### =============== ### =============== ###
###   DEFINITION OF INITIAL PARAMETERS  ###
### =============== ### =============== ###

# N   = 128
# η   = 1.0
# τ   = 6
# rep = 1
# T   = 8*exp10(-5)
# δ   = 0.35

n   = parse(Int, ARGS[1]) # Number of particles
k   = parse(Float64, ARGS[2]) # non local interactions per particle
Tf  = parse(Int, ARGS[3]) # number of iterations
rep = parse(Int, ARGS[4]) # repetition number

# η    = parse(Float64, ARGS[2]) # Dissipation term
# T    = parse(Float64, ARGS[3]) # temperature, noise

@eval @everywhere N    = $n
@eval @everywhere n_nl = $k

@eval @everywhere χ   = 1.25
@eval @everywhere η   = 1.5
@eval @everywhere J   = 0.8
@eval @everywhere T   = 0.01
@eval @everywhere n_t = 6
@eval @everywhere ρ   = 0.3
@eval @everywhere v0  = 0.1

@everywhere pars = InertialParameters(χ, J, η, v0*sqrt(J/χ), ρ, v0, N, T, n_t, n_nl)

# σ = sqrt((2pars.d) * η * T) # noise std deviation ( square root of variance )
@everywhere σ = sqrt((2pars.d) * η * T * pars.dt) # noise std deviation ( square root of variance )
@everywhere L = cbrt(N / pars.ρ) # 3D

### =============== ### =============== ###
### SET UP SYSTEM AND OUTPUT STRUCTURE  ###
### =============== ### =============== ###

Rij = SharedArray{Float64}(N, N)

pos  = SharedArray{Float64}(3N) # particles positions
vel  = SharedArray{Float64}(3N) # array of particles' velocities
v_t  = SharedArray{Float64}(3N) # short-range topological interactions
v_nl = SharedArray{Float64}(3N) # long-range topological interactions
spin = SharedArray{Float64}(3N) # long-range topological interactions

for i in 1:length(pos)
    pos[i]  = 2*rand()*L - L
    vel[i]  = 2*rand() - 1
    spin[i] = 2*rand() - 1
end

for i in 1:3:length(vel)
    norm = sqrt(vel[i]^2 + vel[i+1]^2 + vel[i+2]^2)
    vel[i]   /= norm
    vel[i+1] /= norm
    vel[i+2] /= norm

    vel[i]   *= v0
    vel[i+1] *= v0
    vel[i+2] *= v0

    norm = sqrt(spin[i]^2 + spin[i+1]^2 + spin[i+2]^2)
    spin[i]   /= norm
    spin[i+1] /= norm
    spin[i+2] /= norm

    temp_spin = normalize(cross([spin[i],spin[i+1],spin[i+2]],[vel[i],vel[i+1],vel[i+2]]))

    spin[i]   = temp_spin[1]
    spin[i+1] = temp_spin[2]
    spin[i+2] = temp_spin[3]
end

output_path = set_output_data_structure_lr("EXTENDED_INERTIAL_SPIN", N, η, T, n_nl)

pos_file  = open(output_path * "/pos_$(rep).dat", "w+")
vel_file  = open(output_path * "/vel_$(rep).dat", "w+")
spin_file = open(output_path * "/spin_$(rep).dat", "w+")

# write initial conditions
println("//////// ", 1)
write(pos_file,  pos)
write(vel_file,  vel)
write(spin_file, spin)

### ============== ### ============== ###
###          SYSTEM EVOLUTION         ###
### ============== ### ============== ###

times = [convert(Int, exp10(i)) for i in 0:Tf]

for i in 1:(length(times) - 1)

    for t in (times[i]+1):times[i+1]

        evolve_extended_system_sh(pos, vel, v_t, v_nl, spin, Rij, pars, σ)

        if t % times[i] == 0 || t % div(times[i], exp10(1)) == 0
            println("//////// ", t)
            write(pos_file,  pos)
            write(vel_file,  vel)
            write(spin_file, spin)
        end
    end

end

### ============== ### ============== ### ============== ###

close(pos_file)
close(vel_file)
close(spin_file)

rmprocs(workers())

println("Done all")

### ============== ### ============== ### ============== ###
