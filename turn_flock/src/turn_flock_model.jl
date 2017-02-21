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
    χ  ::Real # generalized moment of intertia
    J  ::Real # strength of alignment interaction
    η  ::Real # friction coefficient
    dt ::Real  # Integration step
    ρ  ::Real # Local initial density
    v0 ::Real  # Speed
    N  ::Real  # Number of particles
    T  ::Real  # generalized temperature
    d  ::Real # unknown parameter (spatial dimension)
    n_c ::Real  # number of topological interactions

    # default constructor
    parameters() = new(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10, 5, 0.15, 1.0)

    # full constructor
    parameters(χ, J, η, dt, ρ, v0, N, T, n_c) = new(χ, J, η, dt, ρ, v0, N, T, 3.0, n_c)
end

### ============== ### ============== ### ============== ###
### SPIN UPDATE
### ============== ### ============== ### ============== ###

function spin_update(pos, vel, v_n, spin, pars, σ)

    # matrix of relative distances between particles
    nij = zeros(Float64, pars.N,pars.N)

    # compute nij entries
    for i in 1:pars.N, j = i:pars.N
        nij[i,j] = norm(pos[i] - pos[j])
        nij[j,i] = nij[i,j]
    end

    # compute local topological interaction
    for i in 1:pars.N
        # v_n[i] = mean([vel[i] for i in findin(nij[:,i], sort(nij[:,i])[2:n_c])])
        v_n[i] = sum( [vel[i] for i in findin(nij[:,i], sort(nij[:,i])[2:pars.n_c+1])] )
    end

    map!((x,y,z) -> x + (pars.J/pars.v0^2)*pars.dt*cross(x, y) - (pars.η/pars.χ)*pars.dt*z + sqrt(pars.dt)*cross(x, (σ/pars.v0)*randn(3)),spin, vel, v_n, spin)

end

### ============== ### ============== ### ============== ###
### SYSTEM EVOLUTION
### ============== ### ============== ### ============== ###

function evolve(pos, vel, vn, spin, pars, σ)

    ### VELOCITY UPDATE
    map!((x,y) -> y + (1./pars.χ)*cross(x,y), vel, spin, vel)

    ### SPIN UPDATE
    spin_update(pos, vel, v_n, spin, pars, σ)

    ### POSITION UPDATE
    map!( (x,v) -> x + pars.dt*v, pos, pos, vel )

end

### ============== ### ============== ### ============== ###

N   = parse(Int, ARGS[1])
η   = parse(Float64, ARGS[2])
T   = parse(Int, ARGS[3]) # number of iterations
rep = parse(Int, ARGS[4])

N   = 64
η   = 0.1
T   = 2
rep = 1

J = 0.8
χ = 1.25

v0 = 0.1

times = [convert(Int, exp10(i)) for i in 0:T]

#      parameters(χ, J, η   ,dt             ,ρ  ,v0  ,N   ,T          ,n_c)
# pars = parameters(χ, J, 60.0, 0.1*sqrt(J/χ), 0.3, 0.1, 512, 8*exp10(-5), 6)
# pars = parameters(χ, J, 15.0, 0.1*sqrt(J/χ), 0.3, 0.1, 512, 8*exp10(-5), 6)
pars = parameters(χ, J, η, v0*sqrt(J/χ), 0.3, v0, N, 8*exp10(-5), 6)

σ = sqrt((2pars.d)*pars.η*pars.T) # noise std deviation ( square root of variance )

# L = sqrt(pars.N / pars.ρ) # 2D
L = cqrt(pars.N / pars.ρ) # 3D


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

# local topological interactions
v_n = [zeros(3) for i in 1:pars.N]

# initialize spins as zero vectors
spin = [zeros(3) for i in 1:pars.N]

# array of  particles' spin (initial random directions)
# spin = pars.v0 * [ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:pars.N ]
# # cross product because dot(si, vi) must = 0, then normalize
# map!((x,y) -> normalize(cross(x,y)), spin, spin, vel)

parent_folder_path = "../TFLOCK_DATA"

folder_path = parent_folder_path * "/data_N_$(pars.N)"

reps_path = folder_path * "/data_N_$(pars.N)_eta_$(ARGS[2])"

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

for i in 1:(length(times) - 1)

    if i > 1

        for t in (times[i]+1):times[i+1]

            evolve(pos, vel, v_n, spin, pars, σ)

            if t % times[i] == 0 || t % times[i-1] == 0
                println("//////// ", t)
                write(pos_file, vcat(pos...))
                write(vel_file, vcat(vel...))
                write(spin_file, vcat(spin...))
            end
        end

    else

        for t in (times[i]+1):times[i+1]

            evolve(pos, vel, v_n, spin, pars, σ)

            if t % times[i] == 0
                println("//////// ", t)
                write(pos_file, vcat(pos...))
                write(vel_file, vcat(vel...))
                write(spin_file, vcat(spin...))
            end
        end

    end

end

close(pos_file)
close(vel_file)
close(spin_file)

println("Done !")
