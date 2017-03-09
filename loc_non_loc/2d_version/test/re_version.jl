### ============== ### ============== ### ============== ###
## Numerical simulations of the Discretized version of the
## inertial spin model Reference: Cavagna...
## Martin Zumaya Hernandez
## 14 / 11 / 2016
### ============== ### ============== ### ============== ###

### ============== ### ============== ### ============== ###
### PARAMETERS ##
### ============== ### ============== ### ============== ###

# immutable parameters
#     N  ::Real  # Number of particles
#     κ  ::Real # Connectivity
#     η  ::Real  # generalized temperature
#     ρ  ::Real # Local initial density
#     ω  ::Real # Interactions relative weight
#     v0 ::Real  # Speed
#     r0 ::Real  # Local interaction range
#     dt ::Real  # Integration step
#
#     # default constructor
#     parameters() = new(1024, 0.0, 0.15, 0.3, 1.0, 10.0, 1.0)
#
#     # full constructor
#     parameters(N, κ, η, ρ, ω, v0, r0) = new(N, κ, η, ρ, ω, v0, r0, 1.0)
# end

### ============== ### ============== ### ============== ###
###                 MODEL IMPLEMENTATION
### ============== ### ============== ### ============== ###

function calc_rij(pos, rij, r0)

    # compute nij entries
    for i in 1:size(rij,1), j in (i+1):size(rij,1)

        d = sqrt((pos[i][1] - pos[j][1])^2 + (pos[i][2] - pos[j][2])^2)

        d < r0 && d > zero(Float64) ? rij[j,i] = one(Float64) : rij[j,i] = zero(Float64)

        # sqrt((pos[i][1] - pos[j][1])^2 + (pos[i][2] - pos[j][2])^2) < r0 && sqrt((pos[i][1] - pos[j][1])^2 + (pos[i][2] - pos[j][2])^2) > zero(Float64) ? rij[j,i] = one(Float64) : rij[i, j] = zero(Float64)

        rij[i,j] = rij[j,i]
    end

end

function part_inter_vel(vel, Nij)
    sum(Nij) > zero(Float64) ? sum(map(*, vel, Nij)) / sum(Nij) : sum(map(*, vel, Nij))
end

function calc_interactions(vel, v_r, v_n, rij, Nij)

    for i in 1:size(Nij, 1)
        v_r[i] = part_inter_vel(vel, view(rij, :, i))
        v_n[i] = part_inter_vel(vel, view(Nij, :, i))
    end

end

function rot_move_part(pos, vel, v_r, v_n)

    prop_angle = atan2(vel[2], vel[1])

    loc_angle = 0.0
    non_loc_angle = 0.0

    i_vx = v_r[1]
    i_vy = v_r[2]

    i_vx != 0.0 || i_vy != 0.0 ? loc_angle = atan2(i_vy, i_vx) - prop_angle : loc_angle = 0.0

    i_vx = v_n[1]
    i_vy = v_n[2]

    i_vx != 0.0 || i_vy != 0.0 ? non_loc_angle = atan2(i_vy, i_vx) - prop_angle : non_loc_angle = 0.0

    # total_angle = pars.ω * loc_angle + (1 - pars.ω) * non_loc_angle + pars.η * (2*rand()*pi - pi);

    total_angle = 0.5 * ( loc_angle + non_loc_angle ) + 0.15 * (2 * rand() * pi - pi);

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
###                 SYSTEM EVOLUTION
### ============== ### ============== ### ============== ###

function evolve(pos, vel, v_r, v_n, rij, Nij, r0)

    calc_rij(pos, rij, r0)
    calc_interactions(vel, v_r, v_n, rij,  Nij)
    map(rot_move_part, pos, vel, v_r, v_n)

end

# @time calc_rij(pos, rij, r0)
# @time part_loc_vel(vel, view(rij, :, 1), pars.r0)
# @time part_non_loc_vel(vel, view(Nij, :, 1))
# @time test(vel, v_r, v_n, rij, Nij, pars.r0)
# @time calc_interactions(pos, vel, v_r, rij, v_n, Nij, pars.r0)
# @time map(rot_move_part, pos, vel, v_r, v_n)

### =============== ### =============== ### =============== ###
### DEFINITION OF INITIAL PARAMETERS
### =============== ### =============== ### =============== ###

reps = 5

κ = parse(Int, ARGS[1]) # average non-local interactions

N  = 128 # number of particles
# N  = 512 # number of particles

# η = 0.15 # noise intensity

ρ  = 0.3 # density
L  = sqrt(N / ρ) # size of box

l = 0.1 # Regimen de velocidad
dt = 1.0 # Time integration step

v0 = 1.0 # particle's speed
T  = 10 # integration time steps

r0 = (v0 * dt) / l

p = κ / (N-1)

# pars = parameters( N, κ, η, ρ, 0.5, v0, r0 )

times = [convert(Int, exp10(i)) for i in 0:3]

### ============== ### ============== ### ============== ###
### OUTPUT
### ============== ### ============== ### ============== ###

folder_path = "../DATA/data_N_$(N)_kappa_$(κ)"

try
    mkdir(folder_path)
catch error
    println("Folder already exists, not created")
end

for r in 1:reps

    println("rep: $(r)")

    ### =============== ### =============== ### =============== ###
    ### INITIALIZATION OF PARTICLES INITIAL POSITIONS AND VELOCIDITES
    ### =============== ### =============== ### =============== ###

    # array of random initial particles' postitions
    # pos = [ [2*rand()*L - L, 2*rand()*L - L, 2*rand()*L - L] for i in 1:pars.N ]
    pos = [ [2*rand()*L - L, 2*rand()*L - L] for i in 1:N ]

    # array of particles' velocities
    # vel = v0 * [ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:N ]
    vel = v0 * [ normalize([2*rand() - 1, 2*rand() - 1]) for i in 1:N ]

    # local metric interactions
    # v_r = [zeros(3) for i in 1:N]
    v_r = [zeros(2) for i in 1:N]

    # non local topological interactions
    # v_n = [zeros(3) for i in 1:N]
    v_n = [zeros(2) for i in 1:N]

    rij = zeros(N, N)
    Nij = zeros(N, N)

    # non-local interaction network definition
    for i in 1:N, j in union(1:i-1, i+1:N)
        rand() < p ? Nij[j, i] = 1 : Nij[j, i] = 0
    end

    # println("Ended Init")

    pos_file = open(folder_path * "/pos_$(r).dat", "w+")
    vel_file = open(folder_path * "/vel_$(r).dat", "w+")

    ### ============== ### ============== ### ============== ###
    ### SYSTEM EVOLUTION
    ### ============== ### ============== ### ============== ###

    for i in 1:(length(times) - 1)

        if i > 1

            for t in (times[i]+1):times[i+1]

                evolve(pos, vel, v_r, v_n, rij, Nij, r0)

                if t % times[i] == 0 || t % times[i-1] == 0
                    # println("//////// ", t)
                    write(pos_file, vcat(pos...))
                    write(vel_file, vcat(vel...))
                end
            end

        else

            for t in (times[i]+1):times[i+1]

                evolve(pos, vel, v_r, v_n, rij, Nij, r0)

                if t % times[i] == 0
                    # println("//////// ", t)
                    write(pos_file, vcat(pos...))
                    write(vel_file, vcat(vel...))
                end
            end

        end

    end

    close(pos_file)
    close(vel_file)

    # println("Done rep")

end

println("Done all")
