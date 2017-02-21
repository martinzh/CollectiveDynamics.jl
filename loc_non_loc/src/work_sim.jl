### ============== ### ============== ### ============== ###
## Numerical simulations of the local and non local
## collective motion model in open space
## Martin Zumaya Hernandez
## 21 / 02 / 2017
### ============== ### ============== ### ============== ###

### ============== ### ============== ### ============== ###
###                 MODEL IMPLEMENTATION
### ============== ### ============== ### ============== ###

function calc_rij(pos, r0)

    rij = zeros(Float64, length(pos), length(pos))

    # compute rij entries
    for i in 1:size(rij,1), j in (i+1):size(rij,1)

        # d = sqrt((pos[i][1] - pos[j][1])^2 + (pos[i][2] - pos[j][2])^2)
        d = norm(pos[i] - pos[j])

        d < r0 && d > zero(Float64) ? rij[j,i] = one(Float64) : rij[j,i] = zero(Float64)

        rij[i,j] = rij[j,i]
    end

    return rij
end

### ============== ### ============== ### ============== ###

function set_Nij(p, Nij)

    N = size(Nij, 1)

    for i in 1:N, j in union(1:i-1, i+1:N)
        rand() < p ? Nij[j, i] = one(Float64) : Nij[j, i] = zero(Float64)
    end
end

### ============== ### ============== ### ============== ###

function calc_interactions(vel, v_n, sp_Nij)

    rows = rowvals(sp_Nij)

    for i in 1:size(sp_Nij, 1)

        range = nzrange(sp_Nij, i)

        length(range) > 0 ? v_n[i] = mean([vel[rows[j]] for j in range ] ) : v_n[i] = zeros(2)

    end

end

### ============== ### ============== ### ============== ###

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

    total_angle = 0.5 * (loc_angle + non_loc_angle) + 0.15 * (2 * rand() * pi - pi);

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

function evolve(pos, vel, v_r, v_n, sp_Nij, r0)

    calc_interactions(vel, v_r, sparse(calc_rij(pos, r0)) ) # local
    calc_interactions(vel, v_n, sp_Nij) # non_local

    map(rot_move_part, pos, vel, v_r, v_n)

end

### =============== ### =============== ### =============== ###
### DEFINITION OF INITIAL PARAMETERS
### =============== ### =============== ### =============== ###


N = parse(Int, ARGS[1]) # number of particles
# N  = 1024
# N  = 512

κ = parse(Float64, ARGS[2]) # average non-local interactions
# κ = 2

T = parse(Int, ARGS[3]) # integration time steps

rep = parse(Int, ARGS[4])

# η = 0.15 # noise intensity

ρ  = 0.3 # density
L  = sqrt(N / ρ) # size of box

l = 0.1 # Regimen de velocidad
dt = 1.0 # Time integration step

v0 = 1.0 # particle's speed

r0 = (v0 * dt) / l

p = κ / (N-1)

times = [convert(Int, exp10(i)) for i in 0:T]

### ============== ### ============== ### ============== ###
### OUTPUT
### ============== ### ============== ### ============== ###

parent_folder_path = "../NLOC_DATA"

folder_path = parent_folder_path * "/DATA/data_N_$(N)"

reps_path = folder_path * "/data_N_$(N)_k_$(κ)"

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

# @time evolve(pos, vel, v_r, v_n, sp_Nij, r0)
#
# @time calc_rij(pos, rij, r0)
# sp_rij = sparse(rij)
#
# @time calc_interactions(vel, v_r, sp_rij)
# @time calc_interactions(vel, v_n, sp_Nij)
#
# @time pmap(rot_move_part, pos, vel, v_r, v_n; distributed=false)
# @time pmap(rot_move_part, pos, vel, v_r, v_n)

println("rep: $(rep)")

### =============== ### =============== ### =============== ###
### INITIALIZATION OF PARTICLES INITIAL POSITIONS AND VELOCIDITES
### =============== ### =============== ### =============== ###

# array of random initial particles' postitions
pos = [ [2*rand()*L - L, 2*rand()*L - L] for i in 1:N ]

# array of particles' velocities
vel = v0 * [ normalize([2*rand() - 1, 2*rand() - 1]) for i in 1:N ]

# local metric interactions
v_r = [zeros(2) for i in 1:N]

# non local topological interactions
v_n = [zeros(2) for i in 1:N]

# non-local interaction network definition
Nij = zeros(Float64, N, N)
set_Nij(p, Nij)

sp_Nij = sparse(Nij)

# println("Ended Init")

pos_file = open(reps_path * "/pos_$(rep).dat", "w+")
vel_file = open(reps_path * "/vel_$(rep).dat", "w+")
net_file = open(reps_path * "/net_$(rep).dat", "w+")

write(net_file, Nij)
close(net_file)

### ============== ### ============== ### ============== ###
### SYSTEM EVOLUTION
### ============== ### ============== ### ============== ###

for i in 1:(length(times) - 1)

    if i > 1

        for t in (times[i]+1):times[i+1]

            evolve(pos, vel, v_r, v_n, sp_Nij, r0)

            if t % times[i] == 0 || t % times[i-1] == 0
                # println("//////// ", t)
                write(pos_file, vcat(pos...))
                write(vel_file, vcat(vel...))
            end
        end

    else

        for t in (times[i]+1):times[i+1]

            evolve(pos, vel, v_r, v_n, sp_Nij, r0)

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

println("Done all")
