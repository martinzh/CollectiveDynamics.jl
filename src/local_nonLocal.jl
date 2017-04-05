### ============== ### ============== ### ============== ###
##    Numerical simulations of the local and non local    ##
##    collective motion model in open space               ##
##    Martin Zumaya Hernandez                             ##
##    21 / 02 / 2017                                      ##
### ============== ### ============== ### ============== ###

### =================== ## ## =========================== ###
##            EXTERNAL DEPENDENCIES IMPORTS                ##
### =================== ## ## =========================== ###

 using Quaternions # Needed for 3D space rotations

### ============== ### ============== ### ============== ###
##                       FLOCK TYPE                       ##
### ============== ### ============== ### ============== ###

"""
    LocNonLocFlock(N, L, v0, p, d)
Inertial Flock type
# Constructor Arguments
* N -> number of particles
* L -> size of box
* v0 -> particles speed
* p -> non-local link probability
* d -> spatial dimension

# Fields
* pos  -> particles positions
* vel  -> particles velocities
* v_r  -> local metric ineraction
* v_n  -> non local ineraction
* Nij  -> ineraction network
* sp_Nij  -> sparse ineraction network
"""
type LocNonLocFlock
    pos ::Array{Array{Float64,1},1}
    vel ::Array{Array{Float64,1},1}
    v_r ::Array{Array{Float64,1},1}
    v_n ::Array{Array{Float64,1},1}
    Nij ::Array{Float64,2}
    sp_Nij ::SparseMatrixCSC{Float64,Int64}

    function LocNonLocFlock(N, L, v0, p, d)

        if d == 2
            # array of random initial particles' postitions within a box of size L
            pos = [ [2*rand()*L - L, 2*rand()*L - L] for i in 1:N ]

            # array of particles' velocities
            vel = v0 * [ normalize([2*rand() - 1, 2*rand() - 1]) for i in 1:N ]

            # local metric interactions
            v_r = [zeros(2) for i in 1:N]

            # non local topological interactions
            v_n = [zeros(2) for i in 1:N]

            # non-local interaction network definition
            Nij = zeros(Float64, N, N)

            set_Nij!(p, Nij)
            sp_Nij = sparse(Nij)

        elseif d == 3
            # array of random initial particles' postitions
            pos = [ [2*rand()*L - L, 2*rand()*L - L, 2*rand()*L - L] for i in 1:N ]

            # array of particles' velocities
            vel = v0 * [ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:N ]

            # local metric interactions
            v_r = [zeros(3) for i in 1:N]

            # non local topological interactions
            v_n = [zeros(3) for i in 1:N]

            # non-local interaction network definition
            Nij = zeros(Float64, N, N)

            set_Nij!(p, Nij)
            sp_Nij = sparse(Nij)
        end

        new(pos, vel, v_r, v_n, Nij, sp_Nij)
    end
end

### ============== ### ============== ### ============== ###
##                  SYSTEM'S PARAMETERS                   ##
### ============== ### ============== ### ============== ###

"""
    LocNonLoclParametersParameters(N, κ, ω, η)
Paramaters for the Local + Non-Local Model of Collective Motion
# Fields
* N ::Int64  # Number of particles
* κ ::Float64 # Non-local connectivity per particle
* ω ::Float64 # Relative interactions weight
* η ::Float64 # Noise Intensity

"""
immutable LocNonLocParameters
    N ::Int64  # Number of particles
    κ ::Float64 # Non-local connectivity per particle
    ω ::Float64 # Relative interactions weight
    η ::Float64 # Noise Intensity
    ρ ::Float64 # Initial local density
    l ::Float64 # Speed regime
    v0 ::Float64 # Particle's speed
    dt ::Float64 # Integration step

    # default constructor
    LocNonLocParameters() = new(128, 6.0, 0.5, 0.15, 0.3, 0.1, 1.0, 1.0)

    # full constructor
    LocNonLocParameters(N, κ, ω, ρ, η) = new(N, κ, ω, η, ρ, 0.1, 1.0, 1.0)
end

### ============== ### ============== ### ============== ###
##      ADJACENCY MATRIX BASED ON METRIC DISTANCES        ##
### ============== ### ============== ### ============== ###
"""
    calc_Rij(pos, r0)
Compute the Adjacency Matrix between particles within a certain metric threshold r0.
Entry Rij == 1 if the link between node i and j exists, Rij == 0 otherwise
"""
function calc_Rij(vec, r0)

    Rij = zeros(Float64, length(vec), length(vec))

    # compute Rij entries
    for i in 1:size(Rij,1), j in (i+1):size(Rij,1)

        d = norm(vec[i] - vec[j])
        d < r0 && d > zero(Float64) ? Rij[j,i] = one(Float64) : Rij[j,i] = zero(Float64)

        Rij[i,j] = Rij[j,i]
    end

    return Rij
end

### ============== ### ============== ### ============== ###
##          SET UP NON-LOCAL INTERACTION NETWORK          ##
### ============== ### ============== ### ============== ###
"""
    set_Nij!(p, Nij)
Set up non-local interaction network between particles, each link is present with probability p. Self interactions are not allowed.
Modifies the supplied matrix Nij.
"""
function set_Nij!(p, Nij)

    N = size(Nij, 1)

    for i in 1:N, j in union(1:i-1, i+1:N)
        rand() < p ? Nij[j, i] = one(Float64) : Nij[j, i] = zero(Float64)
    end
end

### ============== ### ============== ### ============== ###
##                    2D SYSTEM SET UP                    ##
### ============== ### ============== ### ============== ###
"""
    set_up_system_2D!(N, L, v0, p)
Dynamical rules for velocity and position update.
# Arguments
* N -> number of particles
* L -> size of box
* v0 -> particles speed
* p -> non-local interaction probability
"""
function set_up_loc_nonLoc_system_2D!(N, L, v0, p)
    # array of random initial particles' postitions within a box of size L
    pos = [ [2*rand()*L - L, 2*rand()*L - L] for i in 1:N ]

    # array of particles' velocities
    vel = v0 * [ normalize([2*rand() - 1, 2*rand() - 1]) for i in 1:N ]

    # local metric interactions
    v_r = [zeros(2) for i in 1:N]

    # non local topological interactions
    v_n = [zeros(2) for i in 1:N]

    # non-local interaction network definition
    Nij = zeros(Float64, N, N)

    set_Nij!(p, Nij)
    sp_Nij = sparse(Nij)

    return pos, vel, v_r, v_n, Nij, sp_Nij
end

### ============== ### ============== ### ============== ###
##                    3D SYSTEM SET UP                    ##
### ============== ### ============== ### ============== ###
"""
    set_up_system_3D!(N, L, v0, p)
Dynamical rules for velocity and position update.
# Arguments
* N -> number of particles
* L -> size of box
* v0 -> particles speed
* p -> non-local interaction probability
"""
function set_up_loc_nonLoc_system_3D!(N, L, v0, p)
    # array of random initial particles' postitions
    pos = [ [2*rand()*L - L, 2*rand()*L - L, 2*rand()*L - L] for i in 1:N ]

    # array of particles' velocities
    vel = v0 * [ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:N ]

    # local metric interactions
    v_r = [zeros(3) for i in 1:N]

    # non local topological interactions
    v_n = [zeros(3) for i in 1:N]

    # non-local interaction network definition
    Nij = zeros(Float64, N, N)

    set_Nij!(p, Nij)
    sp_Nij = sparse(Nij)

    return pos, vel, v_r, v_n, Nij, sp_Nij
end

### ============== ### ============== ### ============== ###
##                COMPUTE INTERACTION SIGNAL              ##
### ============== ### ============== ### ============== ###
"""
    calc_interactions!(vect, v_in, sp_Mij, d)
Compute and store particles interactions in v_in, given by the adjacency matrix sp_Mij
sp_Mij is a sparse matrix.
The parameter d is the spatial dimension of the system
"""
function calc_interactions!(vect, v_in, sp_Mij, d)

    rows = rowvals(sp_Mij)

    for i in 1:size(sp_Mij, 1)

        range = nzrange(sp_Mij, i)

        length(range) > 0 ? v_in[i] = mean([vect[rows[j]] for j in range ] ) : v_in[i] = zeros(d)

    end

end

### ============== ### ============== ### ============== ###
##            2D PARTICLES' DYNAMICAL RULES               ##
### ============== ### ============== ### ============== ###
"""
    rot_move_part_2D!(pos, vel, v_r, v_n, η, ω)
Dynamical rules for velocity and position update.
# Arguments
* pos -> particles' positions
* vel -> particles' velocities
* v_r -> local interactions signal
* v_n -> non local interaction signal
* η -> noise intensity
* ω -> interactions relative weight
"""
function rot_move_part_2D!(pos, vel, v_r, v_n, η, ω)

    prop_angle = atan2(vel[2], vel[1])

    loc_angle = 0.0
    non_loc_angle = 0.0

    i_vx = v_r[1]
    i_vy = v_r[2]

    i_vx != 0.0 || i_vy != 0.0 ? loc_angle = atan2(i_vy, i_vx) - prop_angle : loc_angle = 0.0

    i_vx = v_n[1]
    i_vy = v_n[2]

    i_vx != 0.0 || i_vy != 0.0 ? non_loc_angle = atan2(i_vy, i_vx) - prop_angle : non_loc_angle = 0.0

    total_angle = ω * loc_angle + (1.0 - ω) * non_loc_angle + 0.15 * (2.0 * rand() * pi - pi);

    c = cos(total_angle)
    s = sin(total_angle)

    vx = vel[1]*c - vel[2]*s;
    vy = vel[1]*s + vel[2]*c;

    vel[1] = vx;
    vel[2] = vy;

    pos[1] += vx;
    pos[2] += vy;

end

### ============== ### ============== ### ============== ###

"""
    rot_move_part_2D_MOD!(pos, vel, v_r, v_n, η, ω)
Dynamical rules for velocity and position update.
Computes weighted average vector instead of using angles
# Arguments
* pos -> particles' positions
* vel -> particles' velocities
* v_r -> local interactions signal
* v_n -> non local interaction signal
* η -> noise intensity
* ω -> interactions relative weight
"""
function rot_move_part_2D_MOD!(pos, vel, v_r, v_n, η, ω)

    prop_angle = atan2(vel[2], vel[1])
    signal_angle = 0.0

    signal = ω * v_r + (1.0 - ω) * v_n

    i_vx = signal[1]
    i_vy = signal[2]

    i_vx != 0.0 || i_vy != 0.0 ? signal_angle = atan2(i_vy, i_vx) - prop_angle : signal_angle = 0.0

    total_angle = signal_angle + 0.15 * (2.0 * rand() * pi - pi);

    c = cos(total_angle)
    s = sin(total_angle)

    vx = vel[1]*c - vel[2]*s;
    vy = vel[1]*s + vel[2]*c;

    vel[1] = vx;
    vel[2] = vy;

    pos[1] += vx;
    pos[2] += vy;

end

### ============== ### ============== ### ============== ###
##            3D PARTICLES' DYNAMICAL RULES               ##
### ============== ### ============== ### ============== ###
"""
    rot_move_part_3D_MOD!(pos, vel, v_r, v_n, η, ω)
Dynamical rules for velocity and position update.
Computes weighted average vector instead of using angles.
# Arguments
* pos -> particles' positions
* vel -> particles' velocities
* v_r -> local interactions signal
* v_n -> non local interaction signal
* η -> noise intensity
* ω -> interactions relative weight
"""
function rot_move_part_3D_MOD!(pos, vel, v_r, v_n, η, ω)

    signal = ω * v_r + (1.0 - ω) * v_n

    signal_angle = acos(dot(vel, signal))

    noise = randn(3)
    noise_angle = acos(dot(normalize(noise), vel))

    q_r = qrotation(cross(vel, signal), signal_angle) * Quaternion(vel)

    u_vel = [q_r.v1, q_r.v2, q_r.v3]

    q_r = qrotation(cross(u_vel, noise), η * noise_angle) * q_r

    u_vel = normalize([q_r.v1, q_r.v2, q_r.v3])

    vel[1] = u_vel[1]
    vel[2] = u_vel[2]
    vel[3] = u_vel[3]

    pos[1] += u_vel[1]
    pos[2] += u_vel[2]
    pos[3] += u_vel[3]

end

### ============== ### ============== ### ============== ###

"""
    rot_move_part_3D!(pos, vel, v_r, v_n, η, ω)
Dynamical rules for velocity and position update.
# Arguments
* pos -> particles' positions
* vel -> particles' velocities
* v_r -> local interactions signal
* v_n -> non local interaction signal
* η -> noise intensity
* ω -> interactions relative weight
"""
function rot_move_part_3D!(pos, vel, v_r, v_n, η, ω)

    loc_angle     = ω * acos(dot(vel, v_r))
    non_loc_angle = (1.0 - ω) * acos(dot(vel, v_n))

    noise = randn(3)
    noise_angle = acos(dot(normalize(noise), vel))

    # q_r = qrotation(cross(vel, v_r), loc_angle + η * noise_angle) * Quaternion(vel)
    #
    # u_vel = [q_r.v1, q_r.v2, q_r.v3]
    #
    # q_r = qrotation(cross(u_vel, v_n), non_loc_angle + η * noise_angle) * q_r

    q_r = qrotation(cross(vel, v_r), loc_angle ) * Quaternion(vel)

    u_vel = [q_r.v1, q_r.v2, q_r.v3]

    q_r = qrotation(cross(u_vel, v_n), non_loc_angle ) * q_r

    u_vel = [q_r.v1, q_r.v2, q_r.v3]

    q_r = qrotation(cross(u_vel, noise), η * noise_angle) * q_r

    u_vel = normalize([q_r.v1, q_r.v2, q_r.v3])

    vel[1] = u_vel[1]
    vel[2] = u_vel[2]
    vel[3] = u_vel[3]

    pos[1] += u_vel[1]
    pos[2] += u_vel[2]
    pos[3] += u_vel[3]

end

### ============== ### ============== ### ============== ###
##                  2D SYSTEM EVOLUTION                   ##
### ============== ### ============== ### ============== ###
"""
    evolve_2D!(pos, vel, v_r, v_n, sp_Nij, r0, η, ω)
Time evolution of the system
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
function evolve_2D!(pos, vel, v_r, v_n, sp_Nij, r0, η, ω)

    calc_interactions!(vel, v_r, sparse(calc_Rij(pos, r0)), 2 ) # local
    calc_interactions!(vel, v_n, sp_Nij, 2) # non_local

    map( (p, v, vr, vn) -> rot_move_part_2D!(p, v, vr, vn, η, ω), pos, vel, v_r, v_n )

end

### ============== ### ============== ### ============== ###

"""
    evolve_2D_MOD!(pos, vel, v_r, v_n, sp_Nij, r0, η, ω)
Time evolution of the system
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
function evolve_2D_MOD!(pos, vel, v_r, v_n, sp_Nij, r0, η, ω)

    calc_interactions!(vel, v_r, sparse(calc_Rij(pos, r0)), 2 ) # local
    calc_interactions!(vel, v_n, sp_Nij, 2) # non_local

    map( (p, v, vr, vn) -> rot_move_part_2D_MOD!(p, v, vr, vn, η, ω), pos, vel, v_r, v_n )

end

### ============== ### ============== ### ============== ###
##                  3D SYSTEM EVOLUTION                   ##
### ============== ### ============== ### ============== ###
"""
    evolve_3D!(pos, vel, v_r, v_n, sp_Nij, r0, η, ω)
Time evolution of the system
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
function evolve_3D!(pos, vel, v_r, v_n, sp_Nij, r0, η, ω)

    calc_interactions!(vel, v_r, sparse(calc_Rij(pos, r0)), 3) # local
    calc_interactions!(vel, v_n, sp_Nij, 3) # non_local

    map( (p, v, vr, vn) -> rot_move_part_3D!(p, v, vr, vn, η, ω), pos, vel, v_r, v_n )

end

### ============== ### ============== ### ============== ###

"""
    evolve_3D_MOD!(pos, vel, v_r, v_n, sp_Nij, r0, η, ω)
Time evolution of the system
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
function evolve_3D_MOD!(pos, vel, v_r, v_n, sp_Nij, r0, η, ω)

    calc_interactions!(vel, v_r, sparse(calc_Rij(pos, r0)), 3) # local
    calc_interactions!(vel, v_n, sp_Nij, 3) # non_local

    map( (p, v, vr, vn) -> rot_move_part_3D!(p, v, vr, vn, η, ω), pos, vel, v_r, v_n )

end

### ============== ### ============== ### ============== ###
##                 OUTPUT DATA STRUCTURE                  ##
### ============== ### ============== ### ============== ###
"""
    set_output_data_structure_2D(N, κ, ω)
Set up folders for output data
# Arguments
* N -> numer of particles
* κ -> non-local average conectivity
* ω -> interactions relative weight
"""
function set_output_data_structure_2D(N, κ, ω)

    parent_folder_path = "$(homedir())/art_DATA/NLOC_DATA"
    folder_path        = parent_folder_path * "/DATA/data_N_$(N)"
    reps_path          = folder_path * "/data_N_$(N)_k_$(κ)_w_$(ω)"

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
        mkdir(reps_path)
    catch error
        println("Parameter folder already exists")
    end

    return reps_path
end

### ============== ### ============== ### ============== ###
"""
    set_output_data_structure_2D_MOD(N, κ, ω)
Set up folders for output data
# Arguments
* N -> numer of particles
* κ -> non-local average conectivity
* ω -> interactions relative weight
"""
function set_output_data_structure_2D_MOD(N, κ, ω)

    parent_folder_path = "$(homedir())/art_DATA/NLOC_DATA_MOD"
    folder_path        = parent_folder_path * "/DATA/data_N_$(N)"
    reps_path          = folder_path * "/data_N_$(N)_k_$(κ)_w_$(ω)"

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
        mkdir(reps_path)
    catch error
        println("Parameter folder already exists")
    end

    return reps_path
end

### ============== ### ============== ### ============== ###

"""
    set_output_data_structure_3D(N, κ, ω)
Set up folders for output data
# Arguments
* N -> numer of particles
* κ -> non-local average conectivity
* ω -> interactions relative weight
"""
function set_output_data_structure_3D(N, κ, ω)

    parent_folder_path = "$(homedir())/art_DATA/NLOC_DATA_3D"
    folder_path        = parent_folder_path * "/DATA/data_N_$(N)"
    reps_path          = folder_path * "/data_N_$(N)_k_$(κ)_w_$(ω)"

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

    return reps_path
end

### ============== ### ============== ### ============== ###
"""
    set_output_data_structure_3D_MOD(N, κ, ω)
Set up folders for output data
# Arguments
* N -> numer of particles
* κ -> non-local average conectivity
* ω -> interactions relative weight
"""
function set_output_data_structure_3D_MOD(N, κ, ω)

    parent_folder_path = "$(homedir())/art_DATA/NLOC_DATA_MOD_3D"
    folder_path        = parent_folder_path * "/DATA/data_N_$(N)"
    reps_path          = folder_path * "/data_N_$(N)_k_$(κ)_w_$(ω)"

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

    return reps_path
end

### ============== ### ============== ### ============== ###
##                   WHOLE TIME EVOLUTION                 ##
### ============== ### ============== ### ============== ###
"""
    full_time_evolution_2D(pos_file, vel_file, times, flock, r0, η, ω)
System time evolution wrapper
"""
function full_time_evolution_2D(pos_file, vel_file, T, flock, r0, η, ω)

    times = [convert(Int, exp10(i)) for i in 0:T]

    for i in 1:(length(times) - 1)

        if i > 1

            for t in (times[i]+1):times[i+1]

                evolve_2D!(flock.pos, flock.vel, flock.v_r, flock.v_n, flock.sp_Nij, r0, η, ω)

                if t % times[i] == 0 || t % times[i-1] == 0
                    println("//////// ", t)
                    write(pos_file, vcat(flock.pos...))
                    write(vel_file, vcat(flock.vel...))
                end
            end

        else

            for t in (times[i]+1):times[i+1]

                evolve_2D!(flock.pos, flock.vel, flock.v_r, flock.v_n, flock.sp_Nij, r0, η, ω)

                if t % times[i] == 0
                    println("//////// ", t)
                    write(pos_file, vcat(flock.pos...))
                    write(vel_file, vcat(flock.vel...))
                end
            end

        end

    end
end

### ============== ### ============== ### ============== ###
"""
    full_time_evolution_2D_MOD(pos_file, vel_file, times, flock, r0, η, ω)
System time evolution wrapper
"""
function full_time_evolution_2D_MOD(pos_file, vel_file, T, flock, r0, η, ω)

    times = [convert(Int, exp10(i)) for i in 0:T]

    for i in 1:(length(times) - 1)

        if i > 1

            for t in (times[i]+1):times[i+1]

                evolve_2D_MOD!(flock.pos, flock.vel, flock.v_r, flock.v_n, flock.sp_Nij, r0, η, ω)

                if t % times[i] == 0 || t % times[i-1] == 0
                    println("//////// ", t)
                    write(pos_file, vcat(flock.pos...))
                    write(vel_file, vcat(flock.vel...))
                end
            end

        else

            for t in (times[i]+1):times[i+1]

                evolve_2D_MOD!(flock.pos, flock.vel, flock.v_r, flock.v_n, flock.sp_Nij, r0, η, ω)

                if t % times[i] == 0
                    println("//////// ", t)
                    write(pos_file, vcat(flock.pos...))
                    write(vel_file, vcat(flock.vel...))
                end
            end

        end

    end
end

### ============== ### ============== ### ============== ###

"""
    full_time_evolution_3D(pos_file, vel_file, times, flock, r0, η, ω)
System time evolution wrapper
"""
function full_time_evolution_3D(pos_file, vel_file, T, flock, r0, η, ω)

    times = [convert(Int, exp10(i)) for i in 0:T]

    for i in 1:(length(times) - 1)

        if i > 1

            for t in (times[i]+1):times[i+1]

                evolve_3D!(flock.pos, flock.vel, flock.v_r, flock.v_n, flock.sp_Nij, r0, η, ω)

                if t % times[i] == 0 || t % times[i-1] == 0
                    println("//////// ", t)
                    write(pos_file, vcat(flock.pos...))
                    write(vel_file, vcat(flock.vel...))
                end
            end

        else

            for t in (times[i]+1):times[i+1]

                evolve_3D!(flock.pos, flock.vel, flock.v_r, flock.v_n, flock.sp_Nij, r0, η, ω)

                if t % times[i] == 0
                    println("//////// ", t)
                    write(pos_file, vcat(flock.pos...))
                    write(vel_file, vcat(flock.vel...))
                end
            end

        end

    end
end

### ============== ### ============== ### ============== ###

"""
    full_time_evolution_3D_MOD(pos_file, vel_file, times, flock, r0, η, ω)
System time evolution wrapper
"""
function full_time_evolution_3D_MOD(pos_file, vel_file, T, flock, r0, η, ω)

    times = [convert(Int, exp10(i)) for i in 0:T]

    for i in 1:(length(times) - 1)

        if i > 1

            for t in (times[i]+1):times[i+1]

                evolve_3D_MOD!(flock.pos, flock.vel, flock.v_r, flock.v_n, flock.sp_Nij, r0, η, ω)

                if t % times[i] == 0 || t % times[i-1] == 0
                    println("//////// ", t)
                    write(pos_file, vcat(flock.pos...))
                    write(vel_file, vcat(flock.vel...))
                end
            end

        else

            for t in (times[i]+1):times[i+1]

                evolve_3D_MOD!(flock.pos, flock.vel, flock.v_r, flock.v_n, flock.sp_Nij, r0, η, ω)

                if t % times[i] == 0
                    println("//////// ", t)
                    write(pos_file, vcat(flock.pos...))
                    write(vel_file, vcat(flock.vel...))
                end
            end

        end

    end
end
