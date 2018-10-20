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

@everywhere χ   = 1.25
@everywhere η   = 1.5
@everywhere J   = 0.8
@everywhere T   = 0.01
@everywhere n_t = 6
@everywhere ρ   = 0.3
@everywhere v0  = 0.1

@everywhere pars = InertialParameters(χ, J, η, v0*sqrt(J/χ), ρ, v0, N, T, n_t, n_nl)

# σ = sqrt((2pars.d) * η * T) # noise std deviation ( square root of variance )
@everywhere σ = sqrt((2pars.d) * η * T * pars.dt) # noise std deviation ( square root of variance )
@everywhere L = cbrt(N / pars.ρ) # 3D

# @everywhere κ_dist = Poisson(n_nl)

### =============== ### =============== ###
### SET UP SYSTEM AND OUTPUT STRUCTURE  ###
### =============== ### =============== ###

Rij = SharedArray{Float64}(N, N)

k_nl = SharedArray{Int64}(N) # particles positions

pos  = SharedArray{Float64}(3N) # particles positions
vel  = SharedArray{Float64}(3N) # array of particles' velocities
v_t  = SharedArray{Float64}(3N) # short-range topological interactions
v_nl = SharedArray{Float64}(3N) # long-range topological interactions
spin = SharedArray{Float64}(3N) # long-range topological interactions

for i in 1:length(pos)
    pos[i]  = 2*rand()*L - L
    vel[i]  = 2*rand() - 1
    spin[i] = 2*rand() - 1
    v_t[i]  = 0.0
    v_nl[i] = 0.0
end

# for i in 1:3:length(vel)
for i in 0:N-1
    k_nl[i+1] = 0
    # norm = sqrt(vel[i]^2 + vel[i+1]^2 + vel[i+2]^2)
    # vel[i]   /= norm
    # vel[i+1] /= norm
    # vel[i+2] /= norm
    p_vel = v0 * normalize([vel[3i+1],vel[3i+2],vel[3i+3]])

    vel[3i+1] = p_vel[1]
    vel[3i+2] = p_vel[2]
    vel[3i+3] = p_vel[3]

    p_spin = normalize([spin[3i+1],spin[3i+2],spin[3i+3]])

    temp_spin = normalize(cross(p_spin, p_vel))
    # temp_spin = normalize(cross([spin[i],spin[i+1],spin[i+2]],[vel[i],vel[i+1],vel[i+2]]))

    spin[3i+1] = temp_spin[1]
    spin[3i+2] = temp_spin[2]
    spin[3i+3] = temp_spin[3]
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

        evolve_extended_system_sh(pos, vel, v_t, v_nl, spin, k_nl, Rij, pars, σ)

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
