### ============== ### ============== ### ============== ###
##  Numerical simulations of the Discretized version of the
##  inertial spin model Reference: Cavagna...
##  Martin Zumaya Hernandez
##  14 / 11 / 2016
### ============== ### ============== ### ============== ###

using Distributions

### ============== ### ============== ### ============== ###
##                      FLOCK TYPES                       ##
### ============== ### ============== ### ============== ###
"""
    InertialFlock(N, L, v0)
Inertial Flock type
# Constructor Arguments
* N -> number of particles
* L -> size of box
* v0 -> particles speed
# Fields
* pos  -> particles positions
* vel  -> particles velocities
* v_t  -> short-range topological ineraction
* Rij  -> relative distances matrix
* spin -> particles spin
"""
type InertialFlock
    pos ::Array{Array{Float64,1},1}
    vel ::Array{Array{Float64,1},1}
    v_t ::Array{Array{Float64,1},1}
    Rij ::Array{Float64,2}
    spin ::Array{Array{Float64,1},1}

    function InertialFlock(N, L, v0)

        Rij = zeros(Float64, N, N)

        # array of random initial particles' postitions
        pos = [ [2*rand()*L - L, 2*rand()*L - L, 2*rand()*L - L] for i in 1:N ]
        # pos = [ [2*rand()*L - L, 2*rand()*L - L, 0.0] for i in 1:N ]

        # array of particles' velocities
        vel = v0 * [ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:N ]
        # vel = v0 * [normalize([1.0, 0.0, 0.0] - [2*δ*rand() - δ, 2*δ*rand() - δ, 2*δ*rand() - δ]) for i in 1:N]

        # short-range topological interactions
        v_t  = [zeros(Float64, 3) for i in 1:N]

        # initialize spins as zero vectors\
        # spin = [zeros(Float64, 3) for i in 1:N]

        # array of  particles' spin (initial random directions)
        spin = [ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:N ]
        # # cross product because dot(si, vi) must = 0, then normalize
        map!((x,y) -> normalize(cross(x,y)), spin, spin, vel)

        new(pos, vel, v_t, Rij, spin)
    end
end

### ============== ### ============== ### ============== ###

"""
    InertialExtFlock(N, L, v0)
Inertial Flock type (with non local interactions)
# Constructor Arguments
* N -> number of particles
* L -> size of box
* v0 -> particles speed
# Fields
* pos  -> particles positions
* vel  -> particles velocities
* v_t  -> short-range topological ineraction
* v_nl  -> long-range topological ineraction
* Rij  -> relative distances matrix
* spin -> particles spin
"""
type InertialExtFlock
    pos ::Array{Array{Float64,1},1}
    vel ::Array{Array{Float64,1},1}
    v_t ::Array{Array{Float64,1},1}
    v_nl ::Array{Array{Float64,1},1}
    Rij ::Array{Float64,2}
    spin ::Array{Array{Float64,1},1}

    function InertialExtFlock(N, L, v0)

        Rij = zeros(Float64, N, N)

        # array of random initial particles' postitions
        pos = [ [2*rand()*L - L, 2*rand()*L - L, 2*rand()*L - L] for i in 1:N ]
        # pos = [ [2*rand()*L - L, 2*rand()*L - L, 0.0] for i in 1:N ]

        # array of particles' velocities
        vel = v0 * [ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:N ]
        # vel = v0 * [normalize([1.0, 0.0, 0.0] - [2*δ*rand() - δ, 2*δ*rand() - δ, 2*δ*rand() - δ]) for i in 1:N]

        # short-range topological interactions
        v_t  = [zeros(Float64, 3) for i in 1:N]

        # long-range topological interactions
        v_nl = [zeros(Float64, 3) for i in 1:N]

        # initialize spins as zero vectors\
        # spin = [zeros(Float64, 3) for i in 1:N]

        # array of  particles' spin (initial random directions)
        spin = [ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:N ]
        # # cross product because dot(si, vi) must = 0, then normalize
        map!((x,y) -> normalize(cross(x,y)), spin, spin, vel)

        new(pos, vel, v_t, v_nl, Rij, spin)
        end
end

### ============== ### ============== ### ============== ###
##                  SYSTEM'S PARAMETERS                   ##
### ============== ### ============== ### ============== ###
"""
    InertialParameters(χ, J, η, dt, ρ, v0, N, T, n_t, n_nl)
Paramaters for the Inertial Spin Model of Collective Motion
# Fields
* χ  ::Float64 # generalized moment of intertia
* J  ::Float64 # strength of alignment interaction
* η  ::Float64 # friction coefficient
* dt ::Float64  # Integration step
* ρ  ::Float64 # Local initial density
* v0 ::Float64  # Speed
* N  ::Int64  # Number of particles
* T  ::Float64  # generalized temperature
* n_t ::Int64  # number of topological interactions
* n_nl ::Float64 # mean of long-range interactions
"""
immutable InertialParameters
    χ  ::Float64 # generalized moment of intertia
    J  ::Float64 # strength of alignment interaction
    η  ::Float64 # friction coefficient
    dt ::Float64  # Integration step
    ρ  ::Float64 # Local initial density
    v0 ::Float64  # Speed
    N  ::Int64  # Number of particles
    T  ::Float64  # generalized temperature
    d  ::Float64 # spatial dimension
    n_t ::Int64  # number of topological interactions
    κ_dist ::Distributions.Poisson{Float64} # out-degree Poisson Distribution
    # n_nl ::Float64 # mean of non-local interactions

    # default constructor
    InertialParameters() = new(1.25, 0.8, 0.3, 1.0, 1.0, 0.1, 10, 0.5, 3.0, 6, Poisson(1.0))

    # full constructor
    InertialParameters(χ, J, η, dt, ρ, v0, N, T, n_t, n_nl) = new(χ, J, η, dt, ρ, v0, N, T, 3.0, n_t, Poisson(n_nl))
end

### ============== ### ============== ### ============== ###
##      RELATIVE DISTACNCES BETWEEN PARTICLES MATRIX      ##
### ============== ### ============== ### ============== ###
"""
    calc_distance_matrix!(pos, Rij)
Compute the relative distances between particles, store in Rij.
"""
function calc_distance_matrix!(vec, Rij)

    N = size(Rij, 1)

    # compute Rij entries
    for i in 1:N, j in (i+1):N
        Rij[i, j] = norm(vec[i] - vec[j])
        Rij[j, i] = Rij[i, j]
    end
end

### ============== ### ============== ### ============== ###
##       COMPUTE LOCAL TOPOLOGICAL INTERACTIONS           ##
### ============== ### ============== ### ============== ###
"""
    calc_topological_interactions!(v_t, Nij, n_t)
Compute and store particles interactions in v_t, given by the adjacency matrix Nij.
#Arguments
* n_t -> number of fixed local topological interactions
"""
function calc_topological_interactions!(vel, v_t, Rij, n_t)

    # compute local topological interaction
    for i in 1:size(Rij,1)
        # v_t[i] = mean( [vel[j] for j in findin(Rij[:,i], sort(Rij[:,i])[2:n_t+1])] )
        v_t[i] = sum( [ vel[j] for j in findin(Rij[:,i], sort(Rij[:,i])[2:n_t+1]) ] )
    end

end

### ============== ### ============== ### ============== ###
##  COMPUTE LOCAL AND NON LOCAL TOPOLOGICAL INTERACTIONS  ##
### ============== ### ============== ### ============== ###
"""
    calc_lr_topological_interactions!(vel, v_t, v_n, Rij, n_t, n_l)
Compute and store particles local topological interactions in v_t and non-local in v_nl, given by the adjacency matrix Nij.
Non-local interactions are chosen at random from the remaining N - (n_t + 1) particles
#Arguments
* n_t -> number of fixed local topological interactions
* n_nl -> number of fixed non-local topological interactions
"""
function calc_lr_topological_interactions!(vel, v_t, v_nl, Rij, n_t, n_nl)

    # compute local topological interaction
    for i in 1:size(Rij,1)
        v_t[i]  = sum( [ vel[j] for j in findin(Rij[:,i], sort(Rij[:,i])[2:n_t+1]) ] )

        # v_nl[i] = sum( [ vel[j] for j in rand(findin(Rij[:,i], sort(Rij[:,i])[n_t+2:end]), n_nl) ] ) # more efficient already knowing n_nl != 0
        n_nl[i] != zero(Int64) ? v_nl[i] = sum( [ vel[j] for j in rand(findin(Rij[:,i], sort(Rij[:,i])[n_t+2:end]), n_nl[i]) ] ) : v_nl[i] = zeros(Float64, length(v_nl[i]))
    end

end

function calc_lr_toplogical_interactions_mean!(vel, v_t, v_nl, Rij, n_t, n_nl)

    # compute local topological interaction
    for i in 1:size(Rij,1)
        v_t[i]  = mean( [ vel[j] for j in findin(Rij[:,i], sort(Rij[:,i])[2:n_t+1]) ] )

        # v_nl[i] = mean( [ vel[j] for j in rand(findin(Rij[:,i], sort(Rij[:,i])[n_t+2:end]), n_nl) ] ) # more efficient already knowing n_nl != 0
        n_nl[i] != zero(Int64) ? v_nl[i] = mean( [ vel[j] for j in rand(findin(Rij[:,i], sort(Rij[:,i])[n_t+2:end]), n_nl[i]) ] ) : v_nl[i] = zeros(Float64, length(v_nl[i]))
    end

end

### ============== ### ============== ### ============== ###
##             DYNAMICAL RULES PER PARTICLE               ##
### ============== ### ============== ### ============== ###
"""
    part_vel_spin_update(vel, v_t, spin, pars, σ)
Dynamical rules for velocity, spin and position update (Updates 1 particle).
# Arguments
* pars -> system's parameters
* σ -> noise related standard deviation
"""
function part_vel_spin_update!(vel, v_t, spin, pars, σ)

    noise = randn(3) * σ

    u_vel = vel + (pars.dt/pars.χ) * cross(spin, vel)

    u_spin =  ( 1.0 - pars.η * pars.dt / pars.χ ) * spin + (pars.J * pars.dt / pars.v0^2) * cross(vel, v_t) + (pars.dt/ pars.v0) * cross(vel, noise)

    spin = u_spin
    vel  = pars.v0  * normalize(u_vel) # codition of constant speed

end

### ============== ### ============== ### ============== ###
##             DYNAMICAL RULES PER PARTICLE               ##
##              WITH NON-LOCAL INTERACTIONS               ##
### ============== ### ============== ### ============== ###
"""
    part_vel_spin_extended_update!(vel, v_t, v_nl, spin, pars, σ)
Dynamical rules for velocity, spin and position update (Updates 1 particle).
# Arguments
* pars -> system's parameters
* σ -> noise related standard deviation
"""
function part_vel_spin_extended_update!(vel, v_t, v_nl, spin, pars, σ)

    noise = randn(3) * σ

    u_vel = vel + (pars.dt/pars.χ) * cross(spin, vel)

    u_spin =  ( 1.0 - pars.η * pars.dt / pars.χ ) * spin[i] + (pars.J * pars.dt / pars.v0^2) * (cross(vel[i], v_t[i]) + cross(vel[i], v_nl[i]) ) + (pars.dt/ pars.v0) * cross(vel[i], noise)

    spin = u_spin
    vel  = pars.v0  * normalize(u_vel) # codition of constant speed

end

### ============== ### ============== ### ============== ###
##             DYNAMICAL RULES (WHOLE SYSTEM)             ##
### ============== ### ============== ### ============== ###
"""
    vel_spin_update!(vel, v_n, spin, pars, σ)
Dynamical rules for velocity, spin and position update.
# Arguments
* pars -> system's parameters
* σ -> noise related standard deviation
"""
function vel_spin_update!(pos, vel, v_t, spin, pars, σ)

    for i in 1:length(vel)
        noise = randn(3) * σ

        u_vel = vel[i] + (pars.dt/pars.χ) * cross(spin[i], vel[i])

        u_spin =  ( 1.0 - pars.η * pars.dt / pars.χ ) * spin[i] + (pars.J * pars.dt / pars.v0^2) * (cross(vel[i], v_t[i])) + (pars.dt/ pars.v0) * cross(vel[i], noise)

        spin[i] = u_spin
        vel[i]  = pars.v0  * normalize(u_vel) # codition of constant speed

        pos[i] += vel[i]*pars.dt
    end
end

### ============== ### ============== ### ============== ###
##             DYNAMICAL RULES (WHOLE SYSTEM)             ##
##              WITH LONG-RANGE INTERACTIONS               ##
### ============== ### ============== ### ============== ###
"""
    vel_spin_nonLocal_update!(vel, v_n, spin, pars, σ)
Dynamical rules for velocity, spin and position update.
# Arguments
* pars -> system's parameters
* σ -> noise related standard deviation
"""
function vel_spin_extended_update!(pos, vel, v_t, v_nl, spin, pars, σ)

    for i in 1:length(vel)
        noise = randn(3) * σ

        u_vel = vel[i] + (pars.dt/pars.χ) * cross(spin[i], vel[i])

        u_spin =  ( 1.0 - pars.η * pars.dt / pars.χ ) * spin[i] + (pars.J * pars.dt / pars.v0^2) * (cross(vel[i], v_t[i]) + cross(vel[i], v_nl[i]) ) + (pars.dt/ pars.v0) * cross(vel[i], noise)

        spin[i] = u_spin
        vel[i]  = pars.v0  * normalize(u_vel) # codition of constant speed

        pos[i] += vel[i]*pars.dt # update positions
    end
end

### ============== ### ============== ### ============== ###
##                     SYSTEM EVOLUTION                   ##
### ============== ### ============== ### ============== ###
"""
    evolve_inertial_system(pos, vel, v_t, spin, Rij, pars, σ)
Time evolution of the system
"""
function evolve_system(pos, vel, v_t, spin, Rij, pars, σ)

    ### COMPUTE RELATIVE DISTANCES
    calc_distance_matrix!(pos, Rij)

    ### COMPUTE INTERACTIONS
    calc_topological_interactions!(vel, v_t, Rij, pars.n_t)

    ### SPIN UPDATE
    vel_spin_update!(pos, vel, v_t, spin, pars, σ)

    # ### POSITION UPDATE
    # map!( (p,v) -> p + pars.dt * v, pos, pos, vel )

end

### ============== ### ============== ### ============== ###
##                     SYSTEM EVOLUTION                   ##
##               WITH NON LOCAL INTERACTIONS              ##
### ============== ### ============== ### ============== ###
"""
    evolve_nonLocal_inertial_system(pos, vel, v_t, v_nl, spin, nij, pars, σ)
Time evolution of the system, with non-local interactions
# Arguments
* pos -> particles' positions
* vel -> particles' velocities
* v_r -> local interactions signal
* v_n -> non local interaction signal
* sp_Nij -> non local sparse adjacency matrix
* r0 -> local metric interaction range
* η -> noise intensity
* ω -> interactions relative weight
"""
function evolve_extended_system(pos, vel, v_t, v_nl, spin, Rij, pars, σ)

    ### COMPUTE RELATIVE DISTANCES
    calc_distance_matrix!(pos, Rij)

    n_nl = rand(pars.κ_dist, pars.N)

    ### COMPUTE INTERACTIONS
    calc_lr_topological_interactions!(vel, v_t, v_nl, Rij, pars.n_t, n_nl)

    ### SPIN UPDATE
    vel_spin_extended_update!(pos, vel, v_t, v_nl, spin, pars, σ)

    # ### POSITION UPDATE
    # map!( (p,v) -> p + pars.dt * v, pos, pos, vel )

end

### ============== ### ============== ### ============== ###
##                 OUTPUT DATA STRUCTURE                  ##
### ============== ### ============== ### ============== ###
"""
    set_output_data_structure(path, N, η, T))
Set up folders for output data
# Arguments
* N -> numer of particles
* η -> friction coefficient
* T -> temperature (noise)
"""
function set_output_data_structure(path, N, η, T)

    parent_folder_path = "$(homedir())/art_DATA/$(path)"
    folder_path        = parent_folder_path * "/DATA/data_N_$(N)"

    reps_path = folder_path * "/eta_$(η)/eta_$(η)_T_$(T)"

    try
        mkdir("$(homedir())/art_DATA")
    catch error
        println("Main data folder already exists")
    end

    try
        mkdir(parent_folder_path)
    catch error
        println("Parent folder already exists")
    end

    try
        mkdir(parent_folder_path * "/DATA")
    catch error
        println("Parent data folder already exists")
    end

    try
        mkdir(folder_path)
    catch error
        println("Folder already exists")
    end

    try
        mkdir(folder_path * "/eta_$(η)")
    catch error
        println("Parent folder already exists")
    end

    try
        mkdir(reps_path)
    catch error
        println("Parameter folder already exists")
    end

    return reps_path
end

### ============== ### ============== ### ============== ###
"""
    set_output_data_structure_inertial_lr(path, N, η, T, n_nl))
Set up folders for output data
# Arguments
* N -> numer of particles
* η -> friction coefficient
* T -> temperature (noise)
* n_nl -> non-local connectivity
"""
function set_output_data_structure_lr(path, N, η, T, n_nl)

    parent_folder_path = "$(homedir())/art_DATA/$(path)"
    folder_path        = parent_folder_path * "/DATA/data_N_$(N)"

    reps_path = folder_path * "/eta_$(η)/eta_$(η)_T_$(T)_nl_$(n_nl)"

    try
        mkdir("$(homedir())/art_DATA")
    catch error
        println("Main data folder already exists")
    end

    try
        mkdir(parent_folder_path)
    catch error
        println("Parent folder already exists")
    end

    try
        mkdir(parent_folder_path * "/DATA")
    catch error
        println("Parent data folder already exists")
    end

    try
        mkdir(folder_path)
    catch error
        println("Folder already exists")
    end

    try
        mkdir(folder_path * "/eta_$(η)")
    catch error
        println("Parent folder already exists")
    end

    try
        mkdir(reps_path)
    catch error
        println("Parameter folder already exists")
    end

    return reps_path
end
