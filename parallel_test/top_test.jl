

### ============== ### ============== ### ============== ###

# addprocs(4)

### ============== ### ============== ### ============== ###

@everywhere using Distributions, Quaternions

# using Distributions
# @everywhere using Quaternions

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

        ri = div(i,3) + 1

        for j in (i+3):3:length(pos)

            rj = div(j,3) + 1

            R_ij[rj,ri] = (pos[i]-pos[j])^2 + (pos[i+1]-pos[j+1])^2 + (pos[i+2]-pos[j+2])^2
        end
    end

end

### ============== ### ============== ### ============== ###
### COMPUTE METRIC SHORT AND LONG RANGE INTERACTIONS
### ============== ### ============== ### ============== ###

@everywhere function compute_topological_interactions(vel::SharedArray,v_r::SharedArray,v_n::SharedArray,R_ij::SharedArray)

    for id in first(localindexes(vel)):3:last(localindexes(vel))

        i = div(id, 3)
        # print(i+1,"|\t")

        v_r[3i+1] = 0.0
        v_r[3i+2] = 0.0
        v_r[3i+3] = 0.0

        v_n[3i+1] = 0.0
        v_n[3i+2] = 0.0
        v_n[3i+3] = 0.0

        # first neighbors
        sh_n = findin(Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N], sort(Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])[2:k_sh+1])

        # next neighbors
        ln_n = findin(Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N], sort(Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])[k_sh+2:end])

        # short-range
        for j in sh_n
            # print(j,"\t")
            v_r[3i+1] += vel[3(j-1)+1] / k_sh
            v_r[3i+2] += vel[3(j-1)+2] / k_sh
            v_r[3i+3] += vel[3(j-1)+3] / k_sh
        end

        # println()
        k_ln = rand(κ_dist)
        # println(k_ln)

        # possible long range
        if k_ln != 0.0
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

        # norm(signal) != zero(Float64) ? signal_angle = acos(dot(p_vel, signal) / norm(signal)) : signal_angle = 0.0

        q_r = Quaternion(zeros(Float64, 3))

        if norm(signal) != zero(Float64)
            signal_angle = acos(dot(p_vel, signal) / norm(signal))
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

    # k_ln = rand(κ_dist, N)

    # @sync begin
    #     for p in workers()
    #         @spawnat p compute_metric_interactions(vel, v_r, v_n, R_ij, k_ln)
    #         @spawnat p update_particles(pos, vel, v_r, v_n)
    #     end
    # end

    @sync begin
        for p in workers()
            @async remotecall_wait(compute_topological_interactions, p, vel, v_r, v_n, R_ij)
            # @async remotecall_wait(update_particles, p, pos, vel, v_r, v_n)
        end
    end

    @sync begin
        for p in workers()
            # @async remotecall_wait(compute_metric_interactions, p, vel, v_r, v_n, R_ij)
            @async remotecall_wait(update_particles, p, pos, vel, v_r, v_n)
        end
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

# @everywhere N = 512
# @everywhere κ = 12
# @everywhere ω = 0.5
#
# rep = 20

@everywhere k_sh = 6

@everywhere ρ = 0.3
@everywhere η = 0.15
@everywhere v0 = 1.0
@everywhere dt = 1.0
@everywhere l = 0.5

@everywhere L  = cbrt(N / ρ) # size of box

@everywhere r0 = ((v0 * dt) / l)^2 # local interaction range

@everywhere κ_dist = Poisson(κ)
# κ_dist = Poisson(κ)

pos = SharedArray{Float64}(3*N) # particles positions
vel = SharedArray{Float64}(3*N) # array of particles' velocities

v_r = SharedArray{Float64}(3*N) # local metric interactions
v_n = SharedArray{Float64}(3*N) # non local topological interactions

R_ij = SharedArray{Float64}(N,N)

### ============== ### ============== ### ============== ###
### SET OUTPUT
### ============== ### ============== ### ============== ###

output_path = set_output_data_structure_lnl("NLOC_P_TOP_3D", N, κ, ω)

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

# write initial conditions
write(pos_file, pos)
write(vel_file, vel)

### ============== ### ============== ### ============== ###

# for t in 1:T
#
#     println("//////// ", t)
#
#     evolve_system(pos, vel, v_r, v_n, R_ij)
#     write(pos_file, pos)
#     write(vel_file, vel)
#
# end
#
# close(pos_file)
# close(vel_file)
#
# println("Done all")

### ============== ### ============== ### ============== ###

times = [convert(Int, exp10(i)) for i in 0:T]

for i in 1:(length(times) - 1)

    if i >= 1

        for t in (times[i]+1):times[i+1]

            evolve_system(pos, vel, v_r, v_n, R_ij)

            if t % times[i] == 0 || t % times[i-1] == 0
                println("//////// ", t)
                write(pos_file, pos)
                write(vel_file, vel)
            end
        end

    else

        for t in (times[i]+1):times[i+1]

            evolve_system(pos, vel, v_r, v_n, R_ij)

            if t % times[i] == 0
                println("//////// ", t)
                write(pos_file, pos)
                write(vel_file, vel)
            end
        end

    end

end

close(pos_file)
close(vel_file)

println("Done all")

### ============== ### ============== ### ============== ###

# i = 5
# sh_n = find(x -> x > zero(Float64), Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])
#
# length(sh_n)
# convert(Float64, length(sh_n))
#
# sum(Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N][sh_n])

# @time calc_Rij(R_ij, pos)
#
# Symmetric(R_ij, :L)
#
# N
# i = 2
# sh_n = find(x -> x > zero(Float64), Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])
# ln_n = find(x -> x < zero(Float64), Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])
#
# # k_ln = rand(κ_dist, N)
#
# # @sync begin
# #     for p in workers()
# #         @spawnat p compute_metric_interactions(vel, v_r, v_n, R_ij, k_ln)
# #         @spawnat p update_particles(pos, vel, v_r, v_n)
# #     end
# # end
#
# @sync begin
#     for p in workers()
#         @async remotecall_wait(compute_metric_interactions, p, vel, v_r, v_n, R_ij)
#         @async remotecall_wait(update_particles, p, pos, vel, v_r, v_n)
#     end
# end
