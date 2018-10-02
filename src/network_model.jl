### ============== ### ============== ### ============== ###
### OUTPUT FOLDER STRUCTURE
### ============== ### ============== ### ============== ###

function set_output_data_structure(path, N, κ, γ)

    parent_folder_path = "$(homedir())/SIM_DATA/$(path)"
    folder_path        = parent_folder_path * "/DATA/data_N_$(N)"
    reps_path          = folder_path * "/data_N_$(N)_k_$(κ)_g_$(γ)"

    try
        mkdir("$(homedir())/SIM_DATA")
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

### ============== ### ============== ###
###   BUILD DIRECTED RANDOM NETWORK   ###
### ============== ### ============== ###

function set_directed_random_network(N, degs)

    neigh = Vector{Vector{Int}}(N)
    Nij   = Vector{Int64}()
    poski = Vector{Int64}()

    neigh[1] = sample(filter(x-> x != 1, collect(1:N)), degs[1], replace=false, ordered = true)
    # println(1, ":", degs[1], " |", neigh[1])

    push!(poski, 1)
    push!(Nij, degs[1])
    for j in neigh[1]
        push!(Nij, j)
    end

    for i in 2:N
    	neigh[i] = sample(filter(x-> x != i, collect(1:N)), degs[i], replace=false, ordered = true)
        k_last = poski[endof(poski)]
        # println(i, ":", degs[i], " |", neigh[i], " |", k_last)
        if degs[i] > 0
            push!(poski, k_last + degs[i-1] + 1)
            push!(Nij, degs[i])
            for j in neigh[i]
                push!(Nij, j)
            end
        else
            push!(Nij, 0)
            push!(poski, k_last + degs[i-1] + 1)
        end
    end

    # println(poski)
    # println(Nij)

    return poski, Nij
end

### ============== ### ============== ###
###     COMPUTE NETWORK INTERACTIONS
### ============== ### ============== ###

function compute_network_interactions_2D(ang::SharedArray, ang_n::SharedArray, Nij::SharedArray, poski::SharedArray)

    for i in first(localindexes(ang)):last(localindexes(ang))

        ang_n[i] = 0.0

        num_phi = 0.0
        den_phi = 0.0

        ki   = Nij[poski[i]]
        init = poski[i]

        for j in 1:ki
            num_phi += cos( ang[ Nij[ init + j] ])
            den_phi += sin( ang[ Nij[ init + j] ])
        end

        ang_n[i] = atan2(den_phi, num_phi)
    end
end

### ============== ### ============== ### ============== ###
### UPDATE PARTICLE'S POSITIONS AND VELOCITIES
### ============== ### ============== ### ============== ###

function update_particles_2D(pos::SharedArray, ang::SharedArray, ang_n::SharedArray, dt::Float64, v0::Float64, γ::Float64, σ::Float64, ξ_dist)

    for i in first(localindexes(ang)):last(localindexes(ang))

        new_angle = ang[i] + γ * dt * sin(ang_n[i] - ang[i]) + σ * rand(ξ_dist)

        n_y = pos[2i] + v0 * dt * sin(new_angle)
        n_x = pos[2i-1] + v0 * dt * cos(new_angle)

        ang[i] = new_angle

        pos[2i] = n_y
        pos[2i-1] = n_x
    end
end

### ============== ### ============== ### ============== ###
### SYSTEM UPDATE (METRIC SHORT-RANGE INTERACTIONS)
### ============== ### ============== ### ============== ###

function evolve_system_2D(pos::SharedArray, ang::SharedArray, ang_n::SharedArray, poski::SharedArray, Nij::SharedArray, dt::Float64, v0::Float64, γ::Float64, σ::Float64, ξ_dist)

    @sync begin
        for p in workers()
            @async remotecall_wait(compute_network_interactions_2D, p, ang, ang_n, Nij, poski)
        end
    end

    @sync begin
        for p in workers()
            @async remotecall_wait(update_particles_2D, p, pos, ang, ang_n, dt, v0, γ, σ, ξ_dist)
        end
    end
end
### ============== ### ============== ### ============== ###
