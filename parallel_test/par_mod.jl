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
### COMPUTE RELATIVE DISTANCES WITH METRIC THRESHOLD
### ============== ### ============== ### ============== ###

function calc_Rij_th(R_ij::SharedArray, pos::SharedArray, r0::Float64)

    @parallel for i in 1:3:length(pos)

        ri = div(i,3) + 1

        for j in (i+3):3:length(pos)

            rj = div(j,3) + 1

            d = (pos[i]-pos[j])^2 + (pos[i+1]-pos[j+1])^2 + (pos[i+2]-pos[j+2])^2
            d <= r0 && d > 0.0 ? R_ij[rj,ri] = 1.0 : R_ij[rj,ri] = -1.0
            # d <= r0 ? R_ij[rj,ri] = 1.0 : R_ij[rj,ri] = -1.0
        end
    end

end

### ============== ### ============== ### ============== ###
### COMPUTE RAW RELATIVE DISTANCES
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

function compute_metric_interactions(vel::SharedArray,v_r::SharedArray,v_n::SharedArray,R_ij::SharedArray)

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
        # k_sh = sum(Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N][sh_n])

        ln_n = find(x -> x < zero(Float64), Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])
        # ln_n = deleteat!(collect(1:N), sh_n)

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

        # if isempty(ln_n) == false && k_ln != 0
        if k_ln != 0.0
            for j in rand(ln_n, k_ln)
                # print(j,"\t")
                v_n[3i+1] += vel[3(j-1)+1] / convert(Float64,k_ln)
                v_n[3i+2] += vel[3(j-1)+2] / convert(Float64,k_ln)
                v_n[3i+3] += vel[3(j-1)+3] / convert(Float64,k_ln)
            end
        end

        # if k_ln[i+1] > 0
        #     for j in rand(ln_n, k_ln[i+1])
        #         # print(j,"\t")
        #         v_n[3i+1] += vel[3(j-1)+1] / k_ln[i+1]
        #         v_n[3i+2] += vel[3(j-1)+2] / k_ln[i+1]
        #         v_n[3i+3] += vel[3(j-1)+3] / k_ln[i+1]
        #     end
        # end

        # println()

    end
end

### ============== ### ============== ### ============== ###
### COMPUTE TOPOLOGICAL SHORT AND LONG RANGE INTERACTIONS
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

function update_particles(pos::SharedArray, vel::SharedArray,v_r::SharedArray,v_n::SharedArray)

    for id in first(localindexes(vel)):3:last(localindexes(vel))

        i = div(id, 3)
        # print(i+1,"|\t")

        signal = ω * normalize([v_r[3i+1] , v_r[3i+2], v_r[3i+3]]) + (1.0 - ω) * normalize([v_n[3i+1] , v_n[3i+2], v_n[3i+3]])
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
### SYSTEM UPDATE (METRIC SHORT RANGE)
### ============== ### ============== ### ============== ###

function evolve_metric_system(pos::SharedArray, vel::SharedArray, v_r::SharedArray, v_n::SharedArray, R_ij::SharedArray, r0::Float64)

    calc_Rij_th(R_ij, pos, r0)

    @sync begin
        for p in workers()
            @async remotecall_wait(compute_metric_interactions, p, vel, v_r, v_n, R_ij)
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
### SYSTEM UPDATE
### ============== ### ============== ### ============== ###

function evolve_topological_system(pos::SharedArray, vel::SharedArray, v_r::SharedArray, v_n::SharedArray, R_ij::SharedArray)

    calc_Rij(R_ij, pos)

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
