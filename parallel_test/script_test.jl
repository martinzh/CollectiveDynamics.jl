

### ============== ### ============== ### ============== ###

# addprocs(4)

### ============== ### ============== ### ============== ###

@everywhere using Distributions, Quaternions

### ============== ### ============== ### ============== ###
### OUTPUT FOLDER STRUCTURE
### ============== ### ============== ### ============== ###

function set_output_data_structure_lnl(path, N, κ, ω)

    parent_folder_path = "$(homedir())/art_DATA/$(path)"
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
### COMPUTE RELATIVE DISTANCES
### ============== ### ============== ### ============== ###

function calc_Rij(R_ij::SharedArray, pos::SharedArray)

    @parallel for i in 1:3:length(pos)
            for j in (i+3):3:length(pos)
                R_ij[div(j,3)+1,div(i,3)+1] = 0.0
                d = (pos[i]-pos[j])^2 + (pos[i+1]-pos[j+1])^2 + (pos[i+2]-pos[j+2])^2
                d <= r0^2 && d > 0.0 ? R_ij[div(j,3)+1,div(i,3)+1] = 1.0 : R_ij[div(j,3)+1,div(i,3)+1] = -1.0
            end
    end

end

### ============== ### ============== ### ============== ###
### COMPUTE METRIC SHORT AND LONG RANGE INTERACTIONS
### ============== ### ============== ### ============== ###

@everywhere function compute_metric_interactions(vel::SharedArray,v_r::SharedArray,v_n::SharedArray,R_ij::SharedArray)

    for id in first(localindexes(vel)):3:last(localindexes(vel))

        i = div(id, 3)
        # print(i+1,"|\t")

        v_r[3i+1] = 0.0
        v_r[3i+2] = 0.0
        v_r[3i+3] = 0.0

        v_n[3i+1] = 0.0
        v_n[3i+2] = 0.0
        v_n[3i+3] = 0.0

        sh_n = find(x -> x > zero(Float64), Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])
        k_sh = length(sh_n)
        # k_sh = sum(sh_n)

        ln_n = find(x -> x < zero(Float64), Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])

        # short-range
        if isempty(sh_n) == false
            for j in sh_n
                # print(j,"\t")
                v_r[3i+1] += vel[3(j-1)+1] / k_sh
                v_r[3i+2] += vel[3(j-1)+2] / k_sh
                v_r[3i+3] += vel[3(j-1)+3] / k_sh
            end
        end

        # println()
        k_ln = rand(κ_dist)
        # println(k_ln)

        # possible long range
        # if k_ln != 0.0
        if isempty(ln_n) == false && k_ln > 0
            for j in rand(ln_n, k_ln)
                # print(j,"\t")
                v_n[3i+1] += vel[3(j-1)+1] / k_ln
                v_n[3i+2] += vel[3(j-1)+2] / k_ln
                v_n[3i+3] += vel[3(j-1)+3] / k_ln
            end
        end

        # println()

    end
end

### ============== ### ============== ### ============== ###
### UPDATE PARTICLE'S POSITIONS AND VELOCITIES
### ============== ### ============== ### ============== ###

@everywhere function update_particles(pos::SharedArray, vel::SharedArray,v_r::SharedArray,v_n::SharedArray)

    for id in first(localindexes(vel)):3:last(localindexes(vel))

        i = div(id, 3)
        # print(i+1,"|\t")

        signal = ω * [v_r[3i+1] , v_r[3i+2], v_r[3i+3]] + (1.0 - ω) * [v_n[3i+1] , v_n[3i+2], v_n[3i+3]]
        signal_angle = 0.0

        p_vel = [vel[3i+1] , vel[3i+2], vel[3i+3]]

        norm(signal) != zero(Float64) ? signal_angle = acos(dot(p_vel, signal) / norm(signal)) : signal_angle = 0.0

        q_r = Quaternion(zeros(Float64, 3))

        if norm(signal) != zero(Float64)
            q_r = qrotation(cross(p_vel, signal), signal_angle + η * (2.0 * rand() * pi - pi)) * Quaternion(p_vel)
        else
            noise = randn(3)
            q_r = qrotation(cross(p_vel, noise), η * acos(dot(normalize(noise), p_vel)) ) * Quaternion(p_vel)
        end

        u_vel = normalize([q_r.v1, q_r.v2, q_r.v3])

        vel[3i+1] = u_vel[1]
        vel[3i+2] = u_vel[2]
        vel[3i+3] = u_vel[3]

        pos[3i+1] += u_vel[1]
        pos[3i+2] += u_vel[2]
        pos[3i+3] += u_vel[3]

    end
end

### ============== ### ============== ### ============== ###
### SYSTEM UPDATE
### ============== ### ============== ### ============== ###

function evolve_system(pos::SharedArray, vel::SharedArray, v_r::SharedArray, v_n::SharedArray, R_ij::SharedArray)

    calc_Rij(R_ij, pos)

    for p in workers()
        # remotecall_fetch(compute_metric_interactions, p, vel, v_r, v_n, R_ij)
        @spawnat p compute_metric_interactions(vel, v_r, v_n, R_ij)
    end

    for p in workers()
        # remotecall_fetch(update_particles, p, pos, vel, v_r, v_n)
        @spawnat p update_particles(pos, vel, v_r, v_n)
    end
end

### ============== ### ============== ### ============== ###
### SYSTEM'S PARAMETERS
### ============== ### ============== ### ============== ###

# N   = parse(Int64, ARGS[1]) # average non-local interactions
# κ   = parse(Float64, ARGS[2]) # average non-local interactions
# ω   = parse(Float64, ARGS[3]) # interactions relative weight

n   = parse(Int64, ARGS[1]) # average non-local interactions
k   = parse(Float64, ARGS[2]) # average non-local interactions
w   = parse(Float64, ARGS[3]) # interactions relative weight

T   = parse(Int, ARGS[4]) # integration time steps
rep = parse(Int, ARGS[5])

### ============== ### ============== ### ============== ###

@eval @everywhere N = $n
@eval @everywhere κ = $k
@eval @everywhere ω = $w

@everywhere ρ = 0.3
@everywhere η = 0.15
@everywhere v0 = 1.0
@everywhere dt = 1.0
@everywhere l = 0.5

@everywhere L  = cbrt(N / ρ) # size of box

@everywhere r0 = (v0 * dt) / l # local interaction range

@everywhere κ_dist = Poisson(κ)

pos = SharedArray{Float64}(3*N) # particles positions
vel = SharedArray{Float64}(3*N) # array of particles' velocities

v_r = SharedArray{Float64}(3*N) # local metric interactions
v_n = SharedArray{Float64}(3*N) # non local topological interactions

R_ij = SharedArray{Float64}(N,N)

### ============== ### ============== ### ============== ###
### SET OUTPUT
### ============== ### ============== ### ============== ###

output_path = set_output_data_structure_lnl("NLOC_P_MET_3D", N, κ, ω)

pos_file = open(output_path * "/pos_$(rep).dat", "w+")
vel_file = open(output_path * "/vel_$(rep).dat", "w+")

### ============== ### ============== ### ============== ###
### INITIALIZATION
### ============== ### ============== ### ============== ###

for i in 1:length(pos)
    pos[i] = 2*rand()*L - L
    vel[i] = 2*rand() - 1
end

for i in 1:3:length(vel)
    norm = sqrt(vel[i]^2 + vel[i+1]^2 + vel[i+2]^2)
    vel[i] /= norm
    vel[i+1] /= norm
    vel[i+2] /= norm
end

### ============== ### ============== ### ============== ###

for t in 1:T

    println("//////// ", t)

    evolve_system(pos, vel, v_r, v_n, R_ij)
    write(pos_file, pos)
    write(vel_file, vel)

end

close(pos_file)
close(vel_file)

println("Done all")

### ============== ### ============== ### ============== ###
