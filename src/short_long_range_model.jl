### ============== ### ============== ### ============== ###
### OUTPUT FOLDER STRUCTURE
### ============== ### ============== ### ============== ###

function set_output_data_structure(path, N, κ, ω)

    parent_folder_path = "$(homedir())/art_DATA/$(path)"
    folder_path = parent_folder_path * "/DATA/data_N_$(N)"
    reps_path = folder_path * "/data_N_$(N)_k_$(κ)_w_$(ω)"

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
        println("Data folder already exists")
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

function set_output_data_structure(path, N, κ, ω, η)

    parent_folder_path = "$(homedir())/art_DATA/$(path)"
    folder_path = parent_folder_path * "/DATA/data_N_$(N)"
    reps_path = folder_path * "/data_N_$(N)_k_$(κ)_w_$(ω)_e_$(η)"

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
        println("Data folder already exists")
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
### COMPUTE RELATIVE (SQUARED) DISTANCES WITH METRIC THRESHOLD
### ============== ### ============== ### ============== ###

function calc_Rij_th(R_ij::SharedArray, pos::SharedArray, r0::Float64)

    @parallel for i in 1:3:length(pos)
        ri = div(i, 3) + 1

        for j in (i + 3):3:length(pos)
            rj = div(j, 3) + 1

            d =
                (pos[i] - pos[j])^2 +
                (pos[i + 1] - pos[j + 1])^2 +
                (pos[i + 2] - pos[j + 2])^2
            d <= r0 ? R_ij[rj, ri] = 1 : R_ij[rj, ri] = -1
        end
    end
end

### ============== ### ============== ### ============== ###

function calc_Rij_th_2D(R_ij::SharedArray, pos::SharedArray, r0::Float64)

    @parallel for i in 1:2:length(pos)
        ri = div(i, 2) + 1

        for j in (i + 2):2:length(pos)
            rj = div(j, 2) + 1

            d = (pos[i] - pos[j])^2 + (pos[i + 1] - pos[j + 1])^2
            d <= r0 ? R_ij[rj, ri] = 1 : R_ij[rj, ri] = -1
        end
    end
end

### ============== ### ============== ### ============== ###
### COMPUTE RAW RELATIVE (SQUARED) DISTANCES
### ============== ### ============== ### ============== ###

function calc_Rij(R_ij::SharedArray, pos::SharedArray)
    @parallel for i in 1:3:length(pos)
        ri = div(i, 3) + 1
        for j in (i + 3):3:length(pos)
            rj = div(j, 3) + 1
            R_ij[rj, ri] =
                (pos[i] - pos[j])^2 +
                (pos[i + 1] - pos[j + 1])^2 +
                (pos[i + 2] - pos[j + 2])^2
        end
    end
end

### ============== ### ============== ### ============== ###

function calc_Rij_2D(R_ij::SharedArray, pos::SharedArray)
    @parallel for i in 1:2:length(pos)
        ri = div(i, 2) + 1
        for j in (i + 2):2:length(pos)
            rj = div(j, 2) + 1
            R_ij[rj, ri] = (pos[i] - pos[j])^2 + (pos[i + 1] - pos[j + 1])^2
        end
    end
end

### ============== ### ============== ### ============== ###
### COMPUTE METRIC SHORT AND LONG RANGE INTERACTIONS
### ============== ### ============== ### ============== ###

function compute_metric_interactions(
    vel::SharedArray,
    v_r::SharedArray,
    v_n::SharedArray,
    R_ij::SharedArray,
    N::Int64,
    κ_dist,
)

    for id in first(localindexes(vel)):3:last(localindexes(vel))
        i = div(id, 3)
        # print(i+1,"|\t")

        v_r[3i + 1] = 0.0
        v_r[3i + 2] = 0.0
        v_r[3i + 3] = 0.0

        v_n[3i + 1] = 0.0
        v_n[3i + 2] = 0.0
        v_n[3i + 3] = 0.0

        sh_n = find(x -> x == 1, Symmetric(R_ij, :L)[((i * N) + 1):((i + 1) * N)])
        ln_n = find(x -> x == -1, Symmetric(R_ij, :L)[((i * N) + 1):((i + 1) * N)])

        k_sh = length(sh_n)

        # short-range
        if k_sh > 0
            for j in sh_n
                # print(j,"\t")
                v_r[3i + 1] += vel[3(j - 1) + 1]
                v_r[3i + 2] += vel[3(j - 1) + 2]
                v_r[3i + 3] += vel[3(j - 1) + 3]
            end
            v_r[3i + 1] /= k_sh
            v_r[3i + 2] /= k_sh
            v_r[3i + 3] /= k_sh
        end

        k_ln = rand(κ_dist)

        if k_ln > 0
            for j in sample(ln_n, k_ln; replace=true)
                v_n[3i + 1] += vel[3(j - 1) + 1]
                v_n[3i + 2] += vel[3(j - 1) + 2]
                v_n[3i + 3] += vel[3(j - 1) + 3]
            end
            v_n[3i + 1] /= k_ln
            v_n[3i + 2] /= k_ln
            v_n[3i + 3] /= k_ln
        end

    end
end

### ============== ### ============== ### ============== ###

function compute_metric_interactions_2D(
    vel::SharedArray,
    v_r::SharedArray,
    v_n::SharedArray,
    R_ij::SharedArray,
    N::Int64,
    κ_dist,
)

    for id in first(localindexes(vel)):2:last(localindexes(vel))
        i = div(id, 2)
        # print(i+1,"|\t")

        v_r[2i + 1] = 0.0
        v_r[2i + 2] = 0.0

        v_n[2i + 1] = 0.0
        v_n[2i + 2] = 0.0

        sh_n = find(x -> x == 1, Symmetric(R_ij, :L)[((i * N) + 1):((i + 1) * N)])
        ln_n = find(x -> x == -1, Symmetric(R_ij, :L)[((i * N) + 1):((i + 1) * N)])

        k_sh = length(sh_n)

        # short-range
        if k_sh > 0
            for j in sh_n
                # print(j,"\t")
                v_r[2i + 1] += vel[2(j - 1) + 1]
                v_r[2i + 2] += vel[2(j - 1) + 2]
            end
            v_r[2i + 1] /= k_sh
            v_r[2i + 2] /= k_sh
        end

        k_ln = rand(κ_dist)

        if k_ln > 0
            for j in sample(ln_n, k_ln; replace=true)
                v_n[2i + 1] += vel[2(j - 1) + 1]
                v_n[2i + 2] += vel[2(j - 1) + 2]
            end
            v_n[2i + 1] /= k_ln
            v_n[2i + 2] /= k_ln
        end

    end
end

### ============== ### ============== ### ============== ###
### COMPUTE TOPOLOGICAL SHORT AND LONG RANGE INTERACTIONS
### ============== ### ============== ### ============== ###

function compute_topological_interactions(
    vel::SharedArray,
    v_r::SharedArray,
    v_n::SharedArray,
    R_ij::SharedArray,
    N::Int64,
    k_sh::Int64,
    κ_dist,
)

    for id in first(localindexes(vel)):3:last(localindexes(vel))
        i = div(id, 3)

        v_r[3i + 1] = 0.0
        v_r[3i + 2] = 0.0
        v_r[3i + 3] = 0.0

        v_n[3i + 1] = 0.0
        v_n[3i + 2] = 0.0
        v_n[3i + 3] = 0.0

        neighbors = sortperm(Symmetric(R_ij, :L)[((i * N) + 1):((i + 1) * N)])

        # short-range
        for j in neighbors[2:(k_sh + 1)]
            v_r[3i + 1] += vel[3(j - 1) + 1]
            v_r[3i + 2] += vel[3(j - 1) + 2]
            v_r[3i + 3] += vel[3(j - 1) + 3]
        end

        v_r[3i + 1] /= k_sh
        v_r[3i + 2] /= k_sh
        v_r[3i + 3] /= k_sh

        k_ln = rand(κ_dist)

        # possible long range
        if k_ln > 0
            for j in sample(neighbors[(k_sh + 2):end], k_ln; replace=true)
                v_n[3i + 1] += vel[3(j - 1) + 1]
                v_n[3i + 2] += vel[3(j - 1) + 2]
                v_n[3i + 3] += vel[3(j - 1) + 3]
            end
            v_n[3i + 1] /= k_ln
            v_n[3i + 2] /= k_ln
            v_n[3i + 3] /= k_ln
        end

    end
end

### ============== ### ============== ### ============== ###

function compute_topological_interactions_2D(
    vel::SharedArray,
    v_r::SharedArray,
    v_n::SharedArray,
    R_ij::SharedArray,
    N::Int64,
    k_sh::Int64,
    κ_dist,
)

    for id in first(localindexes(vel)):3:last(localindexes(vel))
        i = div(id, 2)

        v_r[2i + 1] = 0.0
        v_r[2i + 2] = 0.0

        v_n[2i + 1] = 0.0
        v_n[2i + 2] = 0.0

        neighbors = sortperm(Symmetric(R_ij, :L)[((i * N) + 1):((i + 1) * N)])

        # short-range
        for j in neighbors[2:(k_sh + 1)]
            v_r[2i + 1] += vel[2(j - 1) + 1]
            v_r[2i + 2] += vel[2(j - 1) + 2]
        end

        v_r[2i + 1] /= k_sh
        v_r[2i + 2] /= k_sh

        k_ln = rand(κ_dist)

        # possible long range
        if k_ln > 0
            for j in sample(neighbors[(k_sh + 2):end], k_ln; replace=true)
                v_n[2i + 1] += vel[2(j - 1) + 1]
                v_n[2i + 2] += vel[2(j - 1) + 2]
            end
            v_n[2i + 1] /= k_ln
            v_n[2i + 2] /= k_ln
        end
    end
end

### ============== ### ============== ### ============== ###
### UPDATE PARTICLE'S POSITIONS AND VELOCITIES
### ============== ### ============== ### ============== ###

function update_particles(
    pos::SharedArray,
    vel::SharedArray,
    v_r::SharedArray,
    v_n::SharedArray,
    η::Float64,
    ω::Float64,
)

    q_r = Quaternion(zeros(Float64, 3))

    for id in first(localindexes(vel)):3:last(localindexes(vel))
        i = div(id, 3)

        signal =
            ω .* [v_r[3i + 1], v_r[3i + 2], v_r[3i + 3]] +
            (1.0 - ω) .* [v_n[3i + 1], v_n[3i + 2], v_n[3i + 3]]

        p_vel = [vel[3i + 1], vel[3i + 2], vel[3i + 3]]

        if norm(signal) > zero(Float64)

            signal_angle = dot(p_vel, signal) / (norm(signal) * norm(p_vel))

            signal_angle = ifelse(signal_angle < -1, -1, signal_angle)
            signal_angle = ifelse(signal_angle > 1, 1, signal_angle)

            q_r =
                qrotation(
                    cross(p_vel, signal), acos(signal_angle) + η * (2.0 * rand() * pi - pi)
                ) * Quaternion(p_vel)

        else

            noise = randn(3)

            q_r =
                qrotation(cross(p_vel, noise), η * (2.0 * rand() * pi - pi)) *
                Quaternion(p_vel)

        end

        u_vel = normalize([q_r.v1, q_r.v2, q_r.v3])

        vel[3i + 1] = u_vel[1]
        vel[3i + 2] = u_vel[2]
        vel[3i + 3] = u_vel[3]

        pos[3i + 1] += u_vel[1]
        pos[3i + 2] += u_vel[2]
        pos[3i + 3] += u_vel[3]

    end
end

### ============== ### ============== ### ============== ###

function update_particles_2D(
    pos::SharedArray,
    vel::SharedArray,
    v_r::SharedArray,
    v_n::SharedArray,
    η::Float64,
    ω::Float64,
)

    for id in first(localindexes(vel)):2:last(localindexes(vel))
        i = div(id, 2)

        signal = ω .* [v_r[2i + 1], v_r[2i + 2]] + (1.0 - ω) .* [v_n[2i + 1], v_n[2i + 2]]

        p_vel = [vel[2i + 1], vel[2i + 2]]
        n_vel = zeros(Float64, 2)

        if norm(signal) > zero(Float64)

            signal_angle = dot(p_vel, signal) / (norm(signal) * norm(p_vel))

            signal_angle = ifelse(signal_angle < -1, -1, signal_angle)
            signal_angle = ifelse(signal_angle > 1, 1, signal_angle)

            c = cos(signal_angle + η * (2.0 * rand() * pi - pi))
            s = sin(signal_angle + η * (2.0 * rand() * pi - pi))

            n_vel[1] = p_vel[1] * c - p_vel[2] * s
            n_vel[2] = p_vel[1] * s + p_vel[2] * c

        else

            c = cos(η * (2.0 * rand() * pi - pi))
            s = sin(η * (2.0 * rand() * pi - pi))

            n_vel[1] = p_vel[1] * c - p_vel[2] * s
            n_vel[2] = p_vel[1] * s + p_vel[2] * c

        end

        normalize!(n_vel)

        vel[2i + 1] = n_vel[1]
        vel[2i + 2] = n_vel[2]

        pos[2i + 1] += n_vel[1]
        pos[2i + 2] += n_vel[2]

    end
end

### ============== ### ============== ### ============== ###
### SYSTEM UPDATE (METRIC SHORT-RANGE INTERACTIONS)
### ============== ### ============== ### ============== ###

function evolve_metric_system(
    pos::SharedArray,
    vel::SharedArray,
    v_r::SharedArray,
    v_n::SharedArray,
    R_ij::SharedArray,
    r0::Float64,
    N::Int64,
    η::Float64,
    ω::Float64,
    κ_dist,
)

    calc_Rij_th(R_ij, pos, r0)

    @sync begin
        for p in workers()
            @async remotecall_wait(
                compute_metric_interactions, p, vel, v_r, v_n, R_ij, N, κ_dist
            )
        end
    end

    @sync begin
        for p in workers()
            @async remotecall_wait(update_particles, p, pos, vel, v_r, v_n, η, ω)
        end
    end

end

### ============== ### ============== ### ============== ###

function evolve_metric_system_2D(
    pos::SharedArray,
    vel::SharedArray,
    v_r::SharedArray,
    v_n::SharedArray,
    R_ij::SharedArray,
    r0::Float64,
    N::Int64,
    η::Float64,
    ω::Float64,
    κ_dist,
)

    calc_Rij_th(R_ij, pos, r0)

    @sync begin
        for p in workers()
            @async remotecall_wait(
                compute_metric_interactions_2D, p, vel, v_r, v_n, R_ij, N, κ_dist
            )
        end
    end

    @sync begin
        for p in workers()
            @async remotecall_wait(update_particles_2D, p, pos, vel, v_r, v_n, η, ω)
        end
    end

end

### ============== ### ============== ### ============== ###
### SYSTEM UPDATE (TOPOLOGICAL SHORT-RANGE INTERACTIONS)
### ============== ### ============== ### ============== ###

function evolve_topological_system(
    pos::SharedArray,
    vel::SharedArray,
    v_r::SharedArray,
    v_n::SharedArray,
    R_ij::SharedArray,
    N::Int64,
    k_sh::Int64,
    η::Float64,
    ω::Float64,
    κ_dist,
)

    calc_Rij(R_ij, pos)

    @sync begin
        for p in workers()
            @async remotecall_wait(
                compute_topological_interactions, p, vel, v_r, v_n, R_ij, N, k_sh, κ_dist
            )
        end
    end

    @sync begin
        for p in workers()
            @async remotecall_wait(update_particles, p, pos, vel, v_r, v_n, η, ω)
        end
    end

end

### ============== ### ============== ### ============== ###

function evolve_topological_system_2D(
    pos::SharedArray,
    vel::SharedArray,
    v_r::SharedArray,
    v_n::SharedArray,
    R_ij::SharedArray,
    N::Int64,
    k_sh::Int64,
    η::Float64,
    ω::Float64,
    κ_dist,
)

    calc_Rij(R_ij, pos)

    @sync begin
        for p in workers()
            @async remotecall_wait(
                compute_topological_interactions_2D, p, vel, v_r, v_n, R_ij, N, k_sh, κ_dist
            )
        end
    end

    @sync begin
        for p in workers()
            @async remotecall_wait(update_particles_2D, p, pos, vel, v_r, v_n, η, ω)
        end
    end

end
