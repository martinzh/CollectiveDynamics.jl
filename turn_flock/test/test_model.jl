### ============== ### ============== ### ============== ###
## Numerical simulations of the Discretized version of the
## inertial spin model Reference: Cavagna...
## Martin Zumaya Hernandez
## 14 / 11 / 2016
### ============== ### ============== ### ============== ###

### ============== ### ============== ### ============== ###
### PARAMETERS ##
### ============== ### ============== ### ============== ###

immutable parameters
    χ  ::Float64 # generalized moment of intertia
    J  ::Float64 # strength of alignment interaction
    η  ::Float64 # friction coefficient
    dt ::Float64  # Integration step
    ρ  ::Float64 # Local initial density
    v0 ::Float64  # Speed
    N  ::Int64  # Number of particles
    T  ::Float64  # generalized temperature
    d  ::Float64 # unknown parameter (spatial dimension)
    n_c ::Int64  # number of topological interactions

    # default constructor
    parameters() = new(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10, 5, 3.0, 1.0)

    # full constructor
    parameters(χ, J, η, dt, ρ, v0, N, T, n_c) = new(χ, J, η, dt, ρ, v0, N, T, 3.0, n_c)
end

### ============== ### ============== ### ============== ###
### COMPUTE RELATIVE DISTANCES MATRIX
### ============== ### ============== ### ============== ###

function calc_nij(pos, nij)

    N = size(nij, 1)

    # compute nij entries
    for i in 1:N, j in (i+1):N
        nij[i,j] = norm(pos[i] - pos[j])
        nij[j,i] = nij[i,j]
    end
end

### ============== ### ============== ### ============== ###
### COMPUTE TOPOLOGICAN INTERACTIONS
### ============== ### ============== ### ============== ###

function calc_interactions(v_n, nij, n_c)

    # compute local topological interaction
    for i in 1:size(nij,1)
        # v_n[i] = mean( [vel[j] for j in findin(nij[:,i], sort(nij[:,i])[2:n_c+1])] )
        v_n[i] = sum( [ vel[j] for j in findin(nij[:,i], sort(nij[:,i])[2:n_c+1]) ] )
    end

end

### ============== ### ============== ### ============== ###
### VELOCITY UPDATE
### ============== ### ============== ### ============== ###

function part_vel_update(vel, spin, pars)

    for i in 1:length(vel)
        # u_vel = vel[i] + (pars.dt/pars.χ) * cross(spin[i], vel[i])
        u_vel = (pars.dt/pars.χ) * cross(spin[i], vel[i])

        println(u_vel)
        vel[i] += u_vel
    end

end

### ============== ### ============== ### ============== ###
### SPIN UPDATE
### ============== ### ============== ### ============== ###

function part_spin_update(vel, v_n, spin, pars, σ)

    for i in 1:length(vel)

        noise = normalize(randn(3) * σ)
        # noise = randn(3) * σ

        u_vel = (pars.dt/pars.χ) * cross(spin[i], vel[i])

        # u_spin = vel[i] + (pars.J * pars.dt / pars.v0^2) * cross(vel[i], v_n[i]) - (pars.η * pars.dt/ pars.χ) * spin[i] + (sqrt(pars.dt) / pars.v0) * cross(vel[i], noise)
        u_spin = (pars.J * pars.dt / pars.v0^2) * cross(vel[i], v_n[i]) - (pars.η * pars.dt/ pars.χ) * spin[i] + (sqrt(pars.dt) / pars.v0) * cross(vel[i], noise)

        # println(u_spin)

        # spin[i] = normalize(u_spin)
        spin[i] += u_spin
        vel[i] += u_vel

    end
end

### ============== ### ============== ### ============== ###
### SYSTEM EVOLUTION
### ============== ### ============== ### ============== ###

function evolve(pos, vel, vn, spin, nij, pars, σ)

    ### COMPUTE RELATIVE DISTANCES
    calc_nij(pos, nij)

    ### COMPUTE INTERACTIONS
    calc_interactions(v_n, nij, pars.n_c)

    ### VELOCITY UPDATE
    # map(part_vel_update, vel, spin, pars)
    # part_vel_update(vel, spin, pars)

    ### SPIN UPDATE
    # map(part_spin_update, vel, v_n, spin, pars, σ)
    part_spin_update(vel, v_n, spin, pars, σ)

    ### POSITION UPDATE
    map!( (x,v) -> x + pars.dt*v, pos, pos, vel )

end

# [norm(vel[i]) for i in 1:length(vel)]
# [norm(spin[i]) for i in 1:length(spin)]

### ============== ### ============== ### ============== ###

# N   = parse(Int, ARGS[1])
# η   = parse(Float64, ARGS[2])
# τ   = parse(Int, ARGS[3]) # number of iterations
# rep = parse(Int, ARGS[4])

N   = 64
# η   = 0.1

# explota para η = 32

η   = 15.0
τ   = 4
rep = 1

J = 0.8
χ = 1.25

v0 = 0.1

δ = 0.35

# times = [convert(Int, exp10(i)) for i in 0:τ]

#      parameters(χ, J, η   ,dt             ,ρ  ,v0  ,N   ,T          ,n_c)
# pars = parameters(χ, J, 60.0, 0.1*sqrt(J/χ), 0.3, 0.1, 512, 8*exp10(-5), 6)
# pars = parameters(χ, J, 15.0, 0.1*sqrt(J/χ), 0.3, 0.1, 512, 8*exp10(-5), 6)
pars = parameters(χ, J, η, v0*sqrt(J/χ), 0.3, v0, N, 8*exp10(-5), 6)

σ = sqrt((2pars.d)*pars.η*pars.T) # noise std deviation ( square root of variance )

# L = sqrt(pars.N / pars.ρ) # 2D
L = cbrt(pars.N / pars.ρ) # 3D

# pos = [ 2*rand()*L - L for i in 1:3N ] # array of random initial particles' postitions
# vel = v0 * vcat([ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:N ]...) # array of  particles' velocities
# spin = v0 * vcat([ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:N ]...) # array of  particles' velocities
# vcat([ cross([spin[i], spin[i+1], spin[i+2]], [vel[i], vel[i+1], vel[i+2]]) for i in 1:3:3N ] )

### ============== ### ============== ### ============== ###
### SYSTEM INITIALIZATION
### ============== ### ============== ### ============== ###

# array of random initial particles' postitions
pos = [ [2*rand()*L - L, 2*rand()*L - L, 2*rand()*L - L] for i in 1:pars.N ]
# pos = [ [2*rand()*L - L, 2*rand()*L - L, 0.0] for i in 1:pars.N ]

# array of particles' velocities
vel = pars.v0 * [ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:pars.N ]
# vel = pars.v0 * [ normalize([2*rand() - 1, 2*rand() - 1, 0.0]) for i in 1:pars.N ]
# vel = pars.v0 * [normalize([1.0, 0.0, 0.0] - [2*δ*rand() - δ, 2*δ*rand() - δ, 2*δ*rand() - δ]) for i in 1:N]

# local topological interactions
v_n = [zeros(Float64, 3) for i in 1:pars.N]

# initialize spins as zero vectors\
spin = [zeros(Float64, 3) for i in 1:pars.N]

nij = zeros(Float64, N, N)

# array of  particles' spin (initial random directions)
# spin = pars.v0 * [ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:pars.N ]
# # cross product because dot(si, vi) must = 0, then normalize
# map!((x,y) -> normalize(cross(x,y)), spin, spin, vel)

### ============== ### ============== ### ============== ###
### SET UP OUTPUT DATA STRUCTURE
### ============== ### ============== ### ============== ###

parent_folder_path = "../TFLOCK_DATA"
folder_path = parent_folder_path * "/data_N_$(pars.N)"
# reps_path = folder_path * "/data_N_$(pars.N)_eta_$(ARGS[2])"

reps_path = folder_path * "/data_N_$(pars.N)_eta_$(pars.η)"

try
    mkdir(parent_folder_path)
catch error
    println("Parent folder already exists")
end

try
    mkdir(folder_path)
catch error
    println("Folder already exists")
end

try
    mkdir(reps_path)
catch error
    println("Parameter folder already exists")
end

pos_file   = open(reps_path * "/pos_$(rep).dat", "w")
vel_file   = open(reps_path * "/vel_$(rep).dat", "w")
spin_file  = open(reps_path * "/spin_$(rep).dat", "w")

### ============== ### ============== ### ============== ###
### SYSTEM EVOLUTION
### ============== ### ============== ### ============== ###

# println("pos: ",pos[1]," vel: ", vel[1]," spin: ", spin[1], " v0 = ", norm(vel[1]))
# println("pos: ",pos[1]," vel: ", vel[1]," dot: ", dot(vel[1], spin[1]), " v0 = ", norm(vel[1]))

# for i in 1:100
#     evolve(pos, vel, v_n, spin, nij, pars, σ)
#     # println("pos: ",pos[1]," vel: ", vel[1]," spin: ", spin[1], " v0 = ", norm(vel[1]))
#     println("pos: ",pos[1]," vel: ", vel[1]," dot: ", dot(vel[1], spin[1]), " v0 = ", norm(vel[1]))
#
# end

# for i in 1:(length(times) - 1)
#
#     if i > 1
#
#         for t in (times[i]+1):times[i+1]
#
#             evolve(pos, vel, v_n, spin, nij, pars, σ)
#
#             if t % times[i] == 0 || t % times[i-1] == 0
#                 println("//////// ", t)
#                 write(pos_file, vcat(pos...))
#                 write(vel_file, vcat(vel...))
#                 write(spin_file, vcat(spin...))
#             end
#         end
#
#     else
#
#         for t in (times[i]+1):times[i+1]
#
#             evolve(pos, vel, v_n, spin, nij, pars, σ)
#
#             if t % times[i] == 0
#                 println("//////// ", t)
#                 write(pos_file, vcat(pos...))
#                 write(vel_file, vcat(vel...))
#                 write(spin_file, vcat(spin...))
#             end
#         end
#
#     end
#
# end

for t in 1:convert(Int, exp10(τ))
    evolve(pos, vel, v_n, spin, nij, pars, σ)

    println("//////// ", t)

    write(pos_file, vcat(pos...))
    write(vel_file, vcat(vel...))
    write(spin_file, vcat(spin...))
end

close(pos_file)
close(vel_file)
close(spin_file)

println("Done !")
