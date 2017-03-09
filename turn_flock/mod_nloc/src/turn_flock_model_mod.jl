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
    n_t ::Int64  # number of topological interactions
    n_nl ::Int64 # numer of non-local interactions

    # default constructor
    parameters() = new(1.25, 0.8, 0.3, 1.0, 1.0, 0.1, 10, 0.5, 3.0, 6, 3)

    # full constructor
    parameters(χ, J, η, dt, ρ, v0, N, T, n_t, n_nl) = new(χ, J, η, dt, ρ, v0, N, T, 3.0, n_t, n_nl)
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
### COMPUTE TOPOLOGICAL INTERACTIONS
### ============== ### ============== ### ============== ###

function calc_interactions(v_n, nij, n_t)

    # compute local topological interaction
    for i in 1:size(nij,1)
        # v_n[i] = mean( [vel[j] for j in findin(nij[:,i], sort(nij[:,i])[2:n_t+1])] )
        v_n[i] = sum( [ vel[j] for j in findin(nij[:,i], sort(nij[:,i])[2:n_t+1]) ] )
    end

end

function calc_interactions_mod(v_t, v_nl, nij, n_t, n_nl)

    # compute local topological interaction
    for i in 1:size(nij,1)
        # v_n[i] = mean( [vel[j] for j in findin(nij[:,i], sort(nij[:,i])[2:n_t+1])] )
        v_t[i]  = sum( [ vel[j] for j in findin(nij[:,i], sort(nij[:,i])[2:n_t+1]) ] )
        n_nl != zero(Int) ? v_nl[i] = sum( [ vel[j] for j in rand(findin(nij[:,i], sort(nij[:,i])[n_t+2:end]), n_nl) ] ) : v_nl = zeros(3)
    end

end

### ============== ### ============== ### ============== ###
### SPIN UPDATE
### ============== ### ============== ### ============== ###

function part_vel_spin_update(vel, v_n, spin, pars, σ)

    # noise = normalize(randn(3)) * σ
    noise = randn(3) * σ

    u_vel = vel + (pars.dt/pars.χ) * cross(spin, vel)

    u_spin =  ( 1.0 - pars.η * pars.dt / pars.χ ) * spin + (pars.J * pars.dt / pars.v0^2) * cross(vel, v_n) + (pars.dt/ pars.v0) * cross(vel, noise)

    spin = u_spin
    vel  = pars.v0  * normalize(u_vel) # codition of constant speed

end

function vel_spin_update_mod(vel, v_t, v_nl, spin, pars, σ)

    for i in 1:length(vel)
        # noise = normalize(randn(3)) * σ
        noise = randn(3) * σ

        u_vel = vel[i] + (pars.dt/pars.χ) * cross(spin[i], vel[i])

        u_spin =  ( 1.0 - pars.η * pars.dt / pars.χ ) * spin[i] + (pars.J * pars.dt / pars.v0^2) * (cross(vel[i], v_t[i]) + cross(vel[i], v_nl[i]) ) + (pars.dt/ pars.v0) * cross(vel[i], noise)

        spin[i] = u_spin
        vel[i]  = pars.v0  * normalize(u_vel) # codition of constant speed
    end
end

### ============== ### ============== ### ============== ###
### SYSTEM EVOLUTION
### ============== ### ============== ### ============== ###

function evolve_mod(pos, vel, v_t, v_nl, spin, nij, pars, σ)

    ### COMPUTE RELATIVE DISTANCES
    calc_nij(pos, nij)

    ### COMPUTE INTERACTIONS
    calc_interactions_mod(v_t, v_nl, nij, pars.n_t, pars.n_nl)

    ### SPIN UPDATE
    # map( (v, vn, s) -> part_vel_spin_update(v, vn, s, pars, σ), vel, v_n, spin )
    vel_spin_update_mod(vel, v_t, v_nl, spin, pars, σ)

    ### POSITION UPDATE
    map!( (p,v) -> p + pars.dt * v, pos, pos, vel )

end

### ============== ### ============== ### ============== ###

N    = parse(Int, ARGS[1])
η    = parse(Float64, ARGS[2])
T    = parse(Float64, ARGS[3]) # temperature, noise
n_nl = parse(Int, ARGS[4]) # number of iterations
τ    = parse(Int, ARGS[5]) # number of iterations
rep  = parse(Int, ARGS[6])

# explota para η = 32

# N   = 128
# η   = 1.0
# τ   = 6
# rep = 1
# T   = 8*exp10(-5)

ρ   = 0.3
J   = 0.8
χ   = 1.25
v0  = 0.1
δ   = 0.35
n_t = 6
# n_nl = 3

times = [convert(Int, exp10(i)) for i in 0:τ]

pars = parameters(χ, J, η, v0*sqrt(J/χ), ρ, v0, N, T, n_t, n_nl)

σ = sqrt((2pars.d) * η * T) # noise std deviation ( square root of variance )

# L = sqrt(pars.N / pars.ρ) # 2D
L = cbrt(pars.N / pars.ρ) # 3D

### ============== ### ============== ### ============== ###
### SYSTEM INITIALIZATION
### ============== ### ============== ### ============== ###

nij = zeros(Float64, N, N)

# array of random initial particles' postitions
pos = [ [2*rand()*L - L, 2*rand()*L - L, 2*rand()*L - L] for i in 1:pars.N ]
# pos = [ [2*rand()*L - L, 2*rand()*L - L, 0.0] for i in 1:pars.N ]

# array of particles' velocities
vel = pars.v0 * [ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:pars.N ]
# vel = pars.v0 * [normalize([1.0, 0.0, 0.0] - [2*δ*rand() - δ, 2*δ*rand() - δ, 2*δ*rand() - δ]) for i in 1:N]

# local topological interactions
v_t  = [zeros(Float64, 3) for i in 1:pars.N]
v_nl = [zeros(Float64, 3) for i in 1:pars.N]

# initialize spins as zero vectors\
# spin = [zeros(Float64, 3) for i in 1:pars.N]

# array of  particles' spin (initial random directions)
spin = [ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:pars.N ]
# # cross product because dot(si, vi) must = 0, then normalize
map!((x,y) -> normalize(cross(x,y)), spin, spin, vel)

### ============== ### ============== ### ============== ###
### SET UP OUTPUT DATA STRUCTURE
### ============== ### ============== ### ============== ###

parent_folder_path = "$(homedir())/art_DATA/TFLOCK_DATA_MOD"
folder_path        = parent_folder_path * "/DATA/data_N_$(N)"

reps_path = folder_path * "/eta_$(ARGS[2])/eta_$(pars.η)_T_$(ARGS[3])"

try
    mkdir(parent_folder_path)
catch error
    println("Parent folder already exists")
end

try
    mkdir(parent_folder_path * "/DATA")
catch error
    println("Parent folder already exists")
end

try
    mkdir(folder_path)
catch error
    println("Folder already exists")
end

try
    mkdir(folder_path * "/eta_$(ARGS[2])")
catch error
    println("Parent folder already exists")
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

            evolve_mod(pos, vel, v_t, v_nl, spin, nij, pars, σ)

            if t % times[i] == 0 || t % times[i-1] == 0
                println("//////// ", t)
                write(pos_file, vcat(pos...))
                write(vel_file, vcat(vel...))
                write(spin_file, vcat(spin...))
            end
        end

    else

        for t in (times[i]+1):times[i+1]

            evolve_mod(pos, vel, v_t, v_nl, spin, nij, pars, σ)

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

println("Done")
