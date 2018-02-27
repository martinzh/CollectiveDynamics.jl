### ============== ### ============== ### ============== ###
### SHARED ARRAYS VERSION
### COUZIN COLLECTIVE MOTION MODEL
### ============== ### ============== ### ============== ###

@everywhere using Quaternions

### ============ OUTPUT FOLDER STRUCTURE ============ ###

@everywhere function set_output_data_structure_lnl(path, N, o, a)

    parent_folder_path = joinpath(homedir(),"art_DATA",path)
    folder_path        = joinpath(parent_folder_path,"DATA","data_N_$(N)")
    reps_path          = joinpath(folder_path,"data_N_$(N)_o_$(o)_a_$(a)")

    try
        mkdir(joinpath(homedir(),"art_DATA"))
        println(joinpath(homedir(),"art_DATA"))
    catch error
        println("Main data folder already exists")
    end

    try
        mkdir(parent_folder_path)
        println(parent_folder_path)
    catch error
        println("Parent folder already exists")
    end

    try
        mkdir(joinpath(parent_folder_path,"DATA"))
        println(joinpath(parent_folder_path,"DATA"))
    catch error
        println("Parent folder already exists")
    end

    try
        mkdir(folder_path)
        println(folder_path)
    catch error
        println("Folder already exists")
    end

    try
        mkdir(reps_path)
        println(reps_path)
    catch error
        println("Parameter folder already exists")
    end

    return reps_path
end

### ============ COMPUTE RELATIVE DISTANCES ============ ###

function calc_Rij(Rij::SharedArray, pos::SharedArray)

    N = size(Rij, 1)

    # @parallel for j in 1:N
    #     for i in 1+j:N
    #     Rij[i,j] = norm(pos[i] - pos[j])
    # end

    @parallel for j in 1:3:length(pos)

        rj = div(j,3) + 1

        for i in (j+3):3:length(pos)

            ri = div(i,3) + 1

            R_ij[ri,rj] = sqrt((pos[i]-pos[j])^2 + (pos[i+1]-pos[j+1])^2 + (pos[i+2]-pos[j+2])^2)
        end
    end

end

### ============ COMPUTE INTERACTIONS ============ ###

@everywhere function compute_interactions(v_r::SharedArray, pos::SharedArray, vel::SharedArray, Rij::SharedArray, zor::Float64, zoo::Float64, zoa::Float64)

    F_Rij = Symmetric(Rij, :L)

    for id in first(localindexes(pos)):3:last(localindexes(pos))

        i = div(id,3)
        # println(i)

        v_r[3i+1] = 0.0
        v_r[3i+2] = 0.0
        v_r[3i+3] = 0.0

        # repel_neighbors = find( x-> x > 0.0 && x <= zor, F_Rij[:, i])
        repel_neighbors = find( x-> x > 0.0 && x <= zor, F_Rij[(i*N)+1:(i+1)*N])

        if length(repel_neighbors) > 0

            for j in repel_neighbors
                # v_r[3i+1] -= pos[3(j-1)+1] - pos[3i+1] / F_Rij[j,i]
                # v_r[3i+2] -= pos[3(j-1)+2] - pos[3i+2] / F_Rij[j,i]
                # v_r[3i+3] -= pos[3(j-1)+3] - pos[3i+3] / F_Rij[j,i]
                v_r[3i+1] -= pos[3(j-1)+1] - pos[3i+1] / F_Rij[(i*N) + j]
                v_r[3i+2] -= pos[3(j-1)+2] - pos[3i+2] / F_Rij[(i*N) + j]
                v_r[3i+3] -= pos[3(j-1)+3] - pos[3i+3] / F_Rij[(i*N) + j]
            end

        else

            # orient_neighbors = find( x-> x > zor && x < zoo, F_Rij[:, i])
            # atract_neighbors = find( x-> x > zoo && x < zoa, F_Rij[:, i])

            orient_neighbors = find( x-> x > zor && x < zoo, F_Rij[(i*N)+1:(i+1)*N])
            atract_neighbors = find( x-> x > zoo && x < zoa, F_Rij[(i*N)+1:(i+1)*N])

            v_o = zeros(Float64, 3)
            v_a = zeros(Float64, 3)

            for j in find( x-> x > zor && x < zoo, F_Rij[:, i])
                v_o[1] += vel[3(j-1)+1]
                v_o[2] += vel[3(j-1)+2]
                v_o[3] += vel[3(j-1)+3]
            end

            for j in find( x-> x > zoo && x < zoa, F_Rij[:, i])
                # v_a[1] += pos[3(j-1)+1] - pos[3i+1] / F_Rij[j,i]
                # v_a[2] += pos[3(j-1)+2] - pos[3i+2] / F_Rij[j,i]
                # v_a[3] += pos[3(j-1)+3] - pos[3i+3] / F_Rij[j,i]
                v_a[1] += pos[3(j-1)+1] - pos[3i+1] / F_Rij[(i*N) + j]
                v_a[2] += pos[3(j-1)+2] - pos[3i+2] / F_Rij[(i*N) + j]
                v_a[3] += pos[3(j-1)+3] - pos[3i+3] / F_Rij[(i*N) + j]
            end

            if length(v_a) > 0
                v_r[3i+1] = 0.5 * (v_o[1] + v_a[1])
                v_r[3i+2] = 0.5 * (v_o[2] + v_a[2])
                v_r[3i+3] = 0.5 * (v_o[3] + v_a[3])
            else
                v_r[3i+1] = v_o[1]
                v_r[3i+2] = v_o[2]
                v_r[3i+3] = v_o[3]
            end

        end

    end

end

### ============ UPDATE PARTICLE'S POSITIONS AND VELOCITIES ============ ###

@everywhere function update_particles(v_r::SharedArray, pos::SharedArray, vel::SharedArray, η::Float64, θ::Float64, v0::Float64)

    for id in first(localindexes(pos)):last(localindexes(pos))

        i = div(id, 3)
        # println(i)

        q_r = Quaternion(zeros(Float64, 3))

        signal = [v_r[3i+1],v_r[3i+2],v_r[3i+3]]
        p_vel = [vel[3i+1],vel[3i+2],vel[3i+3]]

        if norm(signal) > 0.0

            signal_angle = dot(p_vel, signal) / (norm(signal)*norm(p_vel))

            signal_angle = ifelse( signal_angle < -1, -1, signal_angle)
            signal_angle = ifelse( signal_angle > 1, 1, signal_angle)

            rot_ang = acos(signal_angle) + η * (2.0 * rand() * pi - pi)

            # q_r = qrotation(cross(vel[i], v_r[i]),  rot_ang) * Quaternion(vel[i])

            if rot_ang < θ
                q_r = qrotation(cross(p_vel, signal),  rot_ang) * Quaternion(p_vel)
            else
                # q_r = qrotation(cross(vel[i], v_r[i]), θ + η * (2.0 * rand() * pi - pi)) * Quaternion(vel[i])
                q_r = qrotation(cross(p_vel, signal), θ) * Quaternion(vel[i])
            end

            u_vel = normalize([q_r.v1, q_r.v2, q_r.v3])

            vel[3i+1] = u_vel[1]
            vel[3i+2] = u_vel[2]
            vel[3i+3] = u_vel[3]

            pos[3i+1] += v0 * u_vel[1]
            pos[3i+2] += v0 * u_vel[2]
            pos[3i+3] += v0 * u_vel[3]

        else

            # noise = randn(3)
            #
            # signal_angle = dot(noise, vel[i]) / (norm(noise)*norm(vel[i]))
            #
            # signal_angle = ifelse( signal_angle < -1, -1, signal_angle)
            # signal_angle = ifelse( signal_angle > 1, 1, signal_angle)
            #
            # q_r = qrotation( cross(vel[i], noise), η * acos(signal_angle) ) * Quaternion(vel[i])

            pos[3i+1] += v0 * p_vel[1]
            pos[3i+2] += v0 * p_vel[2]
            pos[3i+3] += v0 * p_vel[3]

        end

        # u_vel = normalize([q_r.v1, q_r.v2, q_r.v3])
        #
        # vel[i] = u_vel
        # pos[i] += v0 * u_vel

    end
end

### ============ SYSTEM UPDATE ============ ###

function evolve_system(pos::SharedArray, vel::SharedArray, v_r::SharedArray, R_ij::SharedArray, zor::Float64, zoo::Float64, zoa::Float64, η::Float64, θ::Float64, v0::Float64)

    calc_Rij(R_ij, pos)

    @sync begin
        for p in workers()
            @async remotecall_wait(compute_interactions, p, v_r, pos, vel, R_ij, zor, zoo, zoa)
        end
    end

    @sync begin
        for p in workers()
            @async remotecall_wait(update_particles, p, v_r, pos, vel, η, θ, v0)
        end
    end

end

### ============ SYSTEM'S PARAMETERS ============ ###

@everywhere dt = 1.0
@everywhere ρ = 0.3
@everywhere l = 0.5
@everywhere v0 = 1.0
@everywhere η = 0.15
@everywhere θ = 40.0

### ============ METRIC BEHAVIORAL THRESHOLDS ============ ###

# zor = ((v0 * dt) / l)^2 # zone of repulsion
# zoo = zor + ((v0 * dt) / l)^2 # zone of orientation
# zoa = zoo + ((v0 * dt) / l)^2 # zone of attraction

n = parse(Int64, ARGS[1])
# N = 512

o = parse(Float64, ARGS[2])
a = parse(Float64, ARGS[3])
# δr = 10.0

T = parse(Int64, ARGS[4])
# T = 1000

rep = parse(Int64, ARGS[5])
# rep = 1

@eval @everywhere N = $n
@eval @everywhere Δo = $o
@eval @everywhere Δa = $a

@everywhere L  = cbrt(N / ρ) # size of box

@everywhere zor = 1.0 # zone of repulsion
@everywhere zoo = zor + Δo # zone of orientation
@everywhere zoa = zoo + Δa # zone of attraction

### ============ SHARED ARRAY INITIALIZATION ============ ###

pos = SharedArray{Float64}(3N) # particles positions
vel = SharedArray{Float64}(3N) # array of particles' velocities
v_r = SharedArray{Float64}(3N) # local metric interactions

Rij = SharedArray{Int64}(N,N)

### ============ RANDOM INITIAL CONDITIONS ============ ###

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

### ============ SET OUTPUT ============ ###

output_path = set_output_data_structure_lnl("COUZIN_3D", N, ARGS[2], ARGS[3])

pos_file = open(joinpath(output_path,"pos_$(rep).dat"), "w+")
vel_file = open(joinpath(output_path,"vel_$(rep).dat"), "w+")

# write initial conditions
println("//////// ", 1)
write(pos_file, pos)
write(vel_file, vel)

### ============ TIME EVOLUTION ============ ###

times = [convert(Int, exp10(i)) for i in 0:T]

for i in 1:(length(times) - 1)

    for t in (times[i]+1):times[i+1]

        # calc_Rij(Rij)
        # compute_interactions(v_r, pos, vel, Symmetric(Rij, :L), zor, zoo, zoa)
        # update_particles(v_r, pos, vel, θ, v0)

        evolve_system(pos, vel, v_r, Rij, zor, zoo, zoa, η, θ, v0)

        if t % times[i] == 0 || t % div(times[i], exp10(1)) == 0
            println("//////// ", t)
            write(pos_file, pos)
            write(vel_file, vel)
        end
    end

end

close(pos_file)
close(vel_file)

rmprocs(workers())

println("Done all")

### ============== ### ============== ### ============== ###

# for t in 1:T
#     println(t)
#     calc_Rij(Rij)
#     compute_interactions(v_r, pos, vel, Symmetric(Rij, :L), zor, zoo, zoa)
#     update_particles(v_r, pos, vel, θ, v0)
#     write(pos_out, vcat(pos...))
#     write(vel_out, vcat(vel...))
# end
