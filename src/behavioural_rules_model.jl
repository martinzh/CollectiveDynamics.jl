### ============== ### ============== ### ============== ###
### SHARED ARRAYS VERSION
### COUZIN COLLECTIVE MOTION MODEL
### ============== ### ============== ### ============== ###

### ============ OUTPUT FOLDER STRUCTURE ============ ###

function set_output_data_structure(path, N, o, a)

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

### ============ COMPUTE RELATIVE (SQUARED) DISTANCES ============ ###

function calc_Rij(Rij::SharedArray, pos::SharedArray)

    @parallel for i in 1:3:length(pos)

        ri = div(i,3) + 1

        for j in (i+3):3:length(pos)

            rj = div(j,3) + 1

            Rij[rj,ri] = sqrt((pos[i]-pos[j])^2 + (pos[i+1]-pos[j+1])^2 + (pos[i+2]-pos[j+2])^2)
        end
    end

end

### ============ COMPUTE INTERACTIONS ============ ###

function compute_interactions(v_int::SharedArray, v_r::SharedArray, v_o::SharedArray, v_a::SharedArray, pos::SharedArray, vel::SharedArray, Rij::SharedArray, N::Int64, zor::Float64, zoo::Float64, zoa::Float64)

    F_Rij = Symmetric(Rij, :L)

    for id in first(localindexes(pos)):3:last(localindexes(pos))

        i = div(id,3)
        i_st = (i*N) # posicion inicial de la columna de distancias de la part i
        # println(i)

        k_int = zeros(Int64, 3) # number of neighbors

        v_int[3i+1] = 0.0
        v_int[3i+2] = 0.0
        v_int[3i+3] = 0.0

        v_r[3i+1] = 0.0
        v_r[3i+2] = 0.0
        v_r[3i+3] = 0.0

        v_o[3i+1] = 0.0
        v_o[3i+2] = 0.0
        v_o[3i+3] = 0.0

        v_a[3i+1] = 0.0
        v_a[3i+2] = 0.0
        v_a[3i+3] = 0.0

        for j in 1:N

            rij = F_Rij[i_st + j]

            if rij > 0.0 && rij <= zor

                v_r[3i+1] -= (pos[3(j-1)+1] - pos[3i+1]) / rij
                v_r[3i+2] -= (pos[3(j-1)+2] - pos[3i+2]) / rij
                v_r[3i+3] -= (pos[3(j-1)+3] - pos[3i+3]) / rij
                k_int[1] += 1

            else
                if rij > zor && rij <= zoo

                    v_o[3i+1] += vel[3(j-1)+1]
                    v_o[3i+2] += vel[3(j-1)+2]
                    v_o[3i+3] += vel[3(j-1)+3]
                    k_int[2] += 1

                elseif rij > zoo && rij <= zoa

                    v_a[3i+1] += (pos[3(j-1)+1] - pos[3i+1]) / rij
                    v_a[3i+2] += (pos[3(j-1)+2] - pos[3i+2]) / rij
                    v_a[3i+3] += (pos[3(j-1)+3] - pos[3i+3]) / rij
                    k_int[3] += 1

                end
            end

        end

        if k_int[1] > 0

            v_int[3i+1] = v_r[3i+1]
            v_int[3i+2] = v_r[3i+2]
            v_int[3i+3] = v_r[3i+3]

            v_norm = sqrt(v_int[3i+1]^2 + v_int[3i+2]^2 + v_int[3i+3]^2)

            v_int[3i+1] /= v_norm
            v_int[3i+2] /= v_norm
            v_int[3i+3] /= v_norm

        else

            if k_int[2] > 0 && k_int[3] > 0

                v_int[3i+1] = 0.5 * (v_o[3i+1] + v_a[3i+1])
                v_int[3i+2] = 0.5 * (v_o[3i+2] + v_a[3i+2])
                v_int[3i+3] = 0.5 * (v_o[3i+3] + v_a[3i+3])

                v_norm = sqrt(v_int[3i+1]^2 + v_int[3i+2]^2 + v_int[3i+3]^2)

                v_int[3i+1] /= v_norm
                v_int[3i+2] /= v_norm
                v_int[3i+3] /= v_norm


            elseif k_int[2] > 0 && k_int[3] == 0

                v_int[3i+1] = v_o[3i+1]
                v_int[3i+2] = v_o[3i+2]
                v_int[3i+3] = v_o[3i+3]

                v_norm = sqrt(v_int[3i+1]^2 + v_int[3i+2]^2 + v_int[3i+3]^2)

                v_int[3i+1] /= v_norm
                v_int[3i+2] /= v_norm
                v_int[3i+3] /= v_norm


            elseif k_int[3] > 0 && k_int[2] == 0

                v_int[3i+1] = v_a[3i+1]
                v_int[3i+2] = v_a[3i+2]
                v_int[3i+3] = v_a[3i+3]

                v_norm = sqrt(v_int[3i+1]^2 + v_int[3i+2]^2 + v_int[3i+3]^2)

                v_int[3i+1] /= v_norm
                v_int[3i+2] /= v_norm
                v_int[3i+3] /= v_norm

            end

        end

    end

end

### ============ UPDATE PARTICLE'S POSITIONS AND VELOCITIES ============ ###

function update_particles(v_r::SharedArray, pos::SharedArray, vel::SharedArray, η::Float64, θ::Float64, v0::Float64)

    q_r = Quaternion(zeros(Float64, 3))

    for id in first(localindexes(pos)):3:last(localindexes(pos))

        i = div(id, 3)
        # println(i)

        signal = [v_r[3i+1],v_r[3i+2],v_r[3i+3]]
        p_vel  = [vel[3i+1],vel[3i+2],vel[3i+3]]

        if norm(signal) > zero(Float64)

            # signal_angle = dot(p_vel, signal) / (norm(signal)*norm(p_vel))
            signal_angle = dot(p_vel, signal)

            signal_angle = ifelse( signal_angle < -1, -1, signal_angle)
            signal_angle = ifelse( signal_angle > 1, 1, signal_angle)

            q_r = qrotation(cross(p_vel, signal),  acos(signal_angle) + η * (2.0 * rand() * pi - pi)) * Quaternion(p_vel)

            # rot_ang = acos(signal_angle) + η * (2.0 * rand() * pi - pi)

            # if rot_ang < θ
            #     q_r = qrotation(cross(p_vel, signal),  rot_ang) * Quaternion(p_vel)
            # else
            #     # q_r = qrotation(cross(vel[i], v_r[i]), θ + η * (2.0 * rand() * pi - pi)) * Quaternion(vel[i])
            #     q_r = qrotation(cross(p_vel, signal), θ) * Quaternion(vel[i])
            # end

        else
            noise = randn(3)
            q_r = qrotation(cross(p_vel, noise),  η * (2.0 * rand() * pi - pi)) * Quaternion(p_vel)
        end

        u_vel = normalize([q_r.v1, q_r.v2, q_r.v3])

        vel[3i+1] = u_vel[1]
        vel[3i+2] = u_vel[2]
        vel[3i+3] = u_vel[3]

        pos[3i+1] += v0 * u_vel[1]
        pos[3i+2] += v0 * u_vel[2]
        pos[3i+3] += v0 * u_vel[3]

    end
end

### ============ SYSTEM UPDATE ============ ###

function evolve_system(pos::SharedArray, vel::SharedArray, v_int::SharedArray, v_r::SharedArray, v_o::SharedArray, v_a::SharedArray, R_ij::SharedArray, N::Int64, zor::Float64, zoo::Float64, zoa::Float64, η::Float64, θ::Float64, v0::Float64)

    calc_Rij(R_ij, pos)

    @sync begin
        for p in workers()
            @async remotecall_wait(compute_interactions, p, v_int, v_r, v_o, v_a, pos, vel, R_ij, N, zor, zoo, zoa)
        end
    end

    @sync begin
        for p in workers()
            @async remotecall_wait(update_particles, p, v_int, pos, vel, η, θ, v0)
        end
    end

end
