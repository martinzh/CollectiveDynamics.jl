### ============== ### ============== ###
### UTILITY FUNCTIONS FOR SVM 3D
### IN PERIODIC BOUNDARY CONDITION
### ============== ### ============== ###

using Quaternions

### ============== ### ============== ###
"""
    Box(L, M)
Simulation Box
# Constructor Arguments
* L -> size of box
* M -> number of cells per dimension

# Fields
* L -> linear size
* M -> number of cells per dimension

* bulk_cells   -> discrete coordinates
* bottom_cells -> discrete coordinates
* top_cells    -> discrete coordinates
* left_wall    -> discrete coordinates
* right_wall   -> discrete coordinates
* front_wall   -> discrete coordinates
* back_wall   -> discrete coordinates
* corners -> discrete coordinates

* center_cell -> real coordinates of cell center
* vel_cell -> average velocity of particles within cell

* p_per_cell -> Dict: cell_id => [particles id within cell]
"""
type Box
    L :: Float64
    M :: Float64

    bulk_cells   :: Array{Array{Int64,1},3}
    bottom_cells :: Array{Array{Int64,1},1}
    top_cells    :: Array{Array{Int64,1},1}
    left_wall    :: Array{Array{Int64,1},1}
    right_wall   :: Array{Array{Int64,1},1}
    front_wall   :: Array{Array{Int64,1},1}
    back_wall :: Array{Array{Int64,1},1}
    corners   :: Array{Array{Int64,1},3}

    p_per_cell :: Dict{Array{Int,1}, Array{Int, 1}}

    vel_cell :: Array{Array{Float64,1},3}
    center_cell :: Array{Array{Float64,1},3}

    function Box(L, M)

        ### ============== ### ============== ###
        ## indices of cells

        bulk_cells = [[i,j,k] for i in 2:M-1, j in 2:M-1, k in 2:M-1]

        bottom_cells = union([[i,j,1] for i in 2:M-1, j in 1:M], [[i,j,1] for i in [1,M], j in 2:M-1])
        top_cells = union([[i,j,M] for i in 2:M-1, j in 1:M], [[i,j,1] for i in [1,M], j in 2:M-1])

        left_wall = union([[1,j,k] for j in 2:M-1, k in 1:M], [[1,j,k] for j in [1,M], k in 2:M-1])
        right_wall = union([[M,j,k] for j in 2:M-1, k in 1:M], [[M,j,k] for j in [1,M], k in 2:M-1])

        front_wall = union([[i,1,k] for i in 2:M-1, k in 1:M], [[i,1,k] for i in [1,M], k in 2:M-1])
        back_wall = union([[i,M,k] for i in 2:M-1, k in 1:M], [[i,M,k] for i in [1,M], k in 2:M-1])

        corners = [[i,j,k] for i in [1,M], j in [1, M], k in [1,M]]

        ### ============== ### ============== ###

        # particles id's in each box
        p_per_cell = Dict{Array{Int,1}, Array{Int, 1}}()

        # average velocity per box
        vel_cell = Array{Array{Float64,1}}(M,M,M)

        # Initialize velocity vector per cell and compute center of each cell
        for i in 1:length(vel_cell)
            vel_cell[i] = zeros(Float64, 3)
        end

        # position of center of each box
        center_cell = Array{Array{Float64,1}}(M,M,M)

        r = linspace(0., L, M+1)
        mp = collect(r[1:length(r) - 1] + 0.5 * step(r))

        for i in 1:M, j in 1:M, k in 1:M
            # println(i,"\t", j, "\t", k)
            center_cell[i,j,k] = [mp[i], mp[j], mp[k]]
        end

        new(L, float(M), bulk_cells, bottom_cells, top_cells, left_wall, right_wall, front_wall, back_wall, corners, p_per_cell, vel_cell, center_cell)
    end
end

### ============== ### ============== ###
"""
    Flock(N, L, v0, p, d)
Inertial Flock type
# Constructor Arguments
* N -> number of particles
* L -> size of box
* dt -> integration time step
* v0 -> particles speed
* r0 -> local interaction range
* η -> noise intensity

# Fields
* N -> number of particles
* v0 -> particles speed
* r0 -> local interaction range
* dt -> integration time step
* η -> noise intensity

* pos -> particles' positions
* vel -> particles' velocities
* v_r -> interaction vector
* k_r -> local connectivity

* p_cell_id -> particle's celll id
"""
type Flock
    N  :: Int64
    v0 :: Float64
    r0 :: Float64
    dt :: Float64
    η :: Float64

    pos :: Array{Array{Float64, 1}}
    vel :: Array{Array{Float64, 1}}
    v_r :: Array{Array{Float64, 1}}
    k_r :: Array{Float64, 1}

    p_cell_id :: Array{Array{Int64,1}}

    function Flock(N, L, dt, v0, r0, η)

        # particles' positions and velocities
        pos = [ [L*rand(), L*rand(), L*rand()] for i in 1:N ]
        vel = v0 * [ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:N ]

        # interaction vector
        v_r = [zeros(3) for i in 1:N]

        # connectivity per particle
        k_r = zeros(Float64, N)

        # box index for each particle in each dimension
        p_cell_id = [zeros(Int,3) for i in 1:N]

        new(N, v0, r0, dt, η, pos, vel, v_r, k_r, p_cell_id)
    end
end

## ============== ### ============== ###
function av_vel_cell(box, vel)
    for k in keys(box.p_per_cell)
        box.vel_cell[k[1], k[2], k[3]] = mean([vel[i] for i in box.p_per_cell[k]])
    end

end

### ============== ### ============== ###
function assign_cell(flock, box, cell_size)

    box.p_per_cell = Dict{Array{Int,1}, Array{Int, 1}}()

    for i in 1:N
        flock.p_cell_id[i] = convert(Array{Int}, div.(floor.(flock.pos[i]), cell_size)) + 1

        haskey(box.p_per_cell, flock.p_cell_id[i]) ? push!(box.p_per_cell[flock.p_cell_id[i]], i) : box.p_per_cell[flock.p_cell_id[i]] = [i]
    end

end

### ============== ### ============== ###
function check_bulk_cells(p_id, cell_id, p_per_cell)

    k_t = Vector{Float64}(27)
    v_t = Vector{Vector{Float64}}(27)

    c = 1

    # println("center cell: ", cell_id)

    for i in -1:1, j in -1:1, k in -1:1
        key = [cell_id[1]+i, cell_id[2]+j, cell_id[3]+k]
        k_t[c], v_t[c] = part_bulk_cell_interaction(p_id, flock, get!(p_per_cell, key, []))
        c += 1
        # println(key, "\t", get!(p_per_cell, key, []))
    end

    sum(k_t) > 0. ? v_t = sum(v_t) ./ sum(k_t) : v_t = zeros(3)

    return v_t
end

### ============== ### ============== ###
function check_corners(p_id, cell_id, p_per_cell, L)

    k_t = Vector{Float64}(27)
    v_t = Vector{Vector{Float64}}(27)

    c = 1
    # println("corner cell: ", cell_id)

    if cell_id[1] == 1 && cell_id[2] == 1 && cell_id[3] == 1
        for i in [M, 1, 2], j in [M, 1, 2], k in [M, 1, 2]
            key = [i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == 1 && cell_id[2] == M && cell_id[3] == 1
        for i in [M, 1, 2], j in [M-1, M, 1], k in [M, 1, 2]
            key = [i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == M && cell_id[2] == 1 && cell_id[3] == 1
        for i in [M-1, M, 1], j in [M, 1, 2], k in [M, 1, 2]
            key = [i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == M && cell_id[2] == M && cell_id[3] == 1
        for i in [M-1, M, 1], j in [M-1, M, 1], k in [M, 1, 2]
            key = [i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == 1 && cell_id[2] == 1 && cell_id[3] == M
        for i in [M, 1, 2], j in [M, 1, 2], k in [M-1, M, 1]
            key = [i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == 1 && cell_id[2] == M && cell_id[3] == M
        for i in [M, 1, 2], j in [M-1, M, 1], k in [M-1, M, 1]
            key = [i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == M && cell_id[2] == 1 && cell_id[3] == M
        for i in [M-1, M, 1], j in [M, 1, 2], k in [M-1, M, 1]
            key = [i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == M && cell_id[2] == M && cell_id[3] == M
        for i in [M-1, M, 1], j in [M-1, M, 1], k in [M-1, M, 1]
            key = [i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    end

    sum(k_t) > 0. ? v_t = sum(v_t) ./ sum(k_t) : v_t = zeros(3)

    return v_t
end

### ============== ### ============== ###
function check_bottom_cells(p_id, cell_id, p_per_cell, L)

    k_t = Vector{Float64}(27)
    v_t = Vector{Vector{Float64}}(27)

    c = 1
    # println("center cell: ", cell_id)

    if in(cell_id[1], 2:M-1) && in(cell_id[2], 2:M-1)
        for i in -1:1, j in -1:1, k in [M, 1, 2]
            key = [cell_id[1]+i, cell_id[2]+j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif in(cell_id[1], 2:M-1) && cell_id[2] == 1
        for i in -1:1, j in [M, 1, 2], k in [M, 1, 2]
            key = [cell_id[1]+i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif in(cell_id[1], 2:M-1) && cell_id[2] == M
        for i in -1:1, j in [M - 1, M, 1], k in [M, 1, 2]
            key = [cell_id[1]+i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == 1 && in(cell_id[2], 2:M-1)
        for i in [M, 1, 2], j in -1:1, k in [M, 1, 2]
            key = [i, cell_id[2] + j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == M && in(cell_id[2], 2:M-1)
        for i in [M - 1, M, 1], j in -1:1, k in [M, 1, 2]
            key = [i, cell_id[2] + j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    end

    sum(k_t) > 0. ? v_t = sum(v_t) ./ sum(k_t) : v_t = zeros(3)

    return v_t
end

### ============== ### ============== ###
function check_top_cells(p_id, cell_id, p_per_cell, L)
    k_t = Vector{Float64}(27)
    v_t = Vector{Vector{Float64}}(27)

    c = 1
    # println("center cell: ", cell_id)

    if in(cell_id[1], 2:M-1) && in(cell_id[2], 2:M-1)
        for i in -1:1, j in -1:1, k in [M-1, M, 1]
            key = [cell_id[1]+i, cell_id[2]+j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif in(cell_id[1], 2:M-1) && cell_id[2] == 1
        for i in -1:1, j in [M, 1, 2], k in [M-1, M, 1]
            key = [cell_id[1]+i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif in(cell_id[1], 2:M-1) && cell_id[2] == M
        for i in -1:1, j in [M - 1, M, 1], k in [M-1, M, 1]
            key = [cell_id[1]+i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == 1 && in(cell_id[2], 2:M-1)
        for i in [M, 1, 2], j in -1:1, k in [M-1, M, 1]
            key = [i, cell_id[2] + j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
        # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == M && in(cell_id[2], 2:M-1)
        for i in [M - 1, M, 1], j in -1:1, k in [M-1, M, 1]
            key = [i, cell_id[2] + j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    end

    sum(k_t) > 0. ? v_t = sum(v_t) ./ sum(k_t) : v_t = zeros(3)

    return v_t
end

### ============== ### ============== ###
function check_left_cells(p_id, cell_id, p_per_cell, L)

    k_t = Vector{Float64}(27)
    v_t = Vector{Vector{Float64}}(27)

    c = 1
    # println("center cell: ", cell_id)

    if in(cell_id[2], 2:M-1) && in(cell_id[3], 2:M-1)
        for i in [M, 1, 2], j in -1:1, k in -1:1
            key = [i, cell_id[2]+j, cell_id[3]+k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
        # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif in(cell_id[2], 2:M-1) && cell_id[3] == 1
        for i in [M, 1, 2], j in -1:1, k in [M, 1, 2]
            key = [i, cell_id[2]+j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif in(cell_id[2], 2:M-1) && cell_id[3] == M
        for i in [M, 1, 2], j in -1:1, k in [M-1, M, 1]
            key = [i, cell_id[2]+i, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[2] == 1 && in(cell_id[3], 2:M-1)
        for i in [M, 1, 2], j in [M, 1, 2], k in -1:1
            key = [i, j, cell_id[3]+k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[2] == M && in(cell_id[3], 2:M-1)
        for i in [M, 1, 2], j in [M-1, M, 1], k in -1:1
            key = [i, j, cell_id[3]+k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    end

    sum(k_t) > 0. ? v_t = sum(v_t) ./ sum(k_t) : v_t = zeros(3)

    return v_t
end

### ============== ### ============== ###
function check_right_cells(p_id, cell_id, p_per_cell, L)

    k_t = Vector{Float64}(27)
    v_t = Vector{Vector{Float64}}(27)

    c = 1
    # println("center cell: ", cell_id)

    if in(cell_id[2], 2:M-1) && in(cell_id[3], 2:M-1)
        for i in [M-1, M, 1], j in -1:1, k in -1:1
            key = [i, cell_id[2]+j, cell_id[3]+k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif in(cell_id[2], 2:M-1) && cell_id[3] == 1
        for i in [M-1, M, 1], j in -1:1, k in [M, 1, 2]
            key = [i, cell_id[2]+j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif in(cell_id[2], 2:M-1) && cell_id[3] == M
        for i in [M-1, M, 1], j in -1:1, k in [M-1, M, 1]
            key = [i, cell_id[2]+j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[2] == 1 && in(cell_id[3], 2:M-1)
        for i in [M-1, M, 1], j in [M, 1, 2], k in -1:1
            key = [i, j, cell_id[3]+k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[2] == M && in(cell_id[3], 2:M-1)
        for i in [M-1, M, 1], j in [M-1, M, 1], k in -1:1
            key = [i, j, cell_id[3]+k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    end

    sum(k_t) > 0. ? v_t = sum(v_t) ./ sum(k_t) : v_t = zeros(3)

    return v_t
end

### ============== ### ============== ###
function check_front_cells(p_id, cell_id, p_per_cell, L)

    k_t = Vector{Float64}(27)
    v_t = Vector{Vector{Float64}}(27)

    c = 1
    # println("center cell: ", cell_id)

    if in(cell_id[1], 2:M-1) && in(cell_id[3], 2:M-1)
        for i in -1:1, j in [M, 1, 2], k in -1:1
            key = [cell_id[1]+i, j, cell_id[3]+k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif in(cell_id[1], 2:M-1) && cell_id[3] == 1
        for i in -1:1, j in [M, 1, 2], k in [M, 1, 2]
            key = [cell_id[1]+i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif in(cell_id[1], 2:M-1) && cell_id[3] == M
        for i in -1:1, j in [M, 1, 2], k in [M-1, M, 1]
            key = [cell_id[1]+i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == 1 && in(cell_id[3], 2:M-1)
        for i in [M, 1, 2], j in [M, 1, 2], k in -1:1
            key = [i, j, cell_id[3]+k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == M && in(cell_id[3], 2:M-1)
        for i in [M-1, M, 1], j in [M, 1, 2], k in -1:1
            key = [i, j, cell_id[3]+k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    end

    sum(k_t) > 0. ? v_t = sum(v_t) ./ sum(k_t) : v_t = zeros(3)

    return v_t
end

### ============== ### ============== ###
function check_back_cells(p_id, cell_id, p_per_cell, L)

    k_t = Vector{Float64}(27)
    v_t = Vector{Vector{Float64}}(27)

    c = 1

    # println("center cell: ", cell_id, " particle id: ", p_id)

    if in(cell_id[1], 2:M-1) && in(cell_id[3], 2:M-1)
        for i in -1:1, j in [M-1, M, 1], k in -1:1
            key = [cell_id[1]+i, j, cell_id[3]+k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif in(cell_id[1], 2:M-1) && cell_id[3] == 1
        for i in -1:1, j in [M-1, M, 1], k in [M, 1, 2]
            key = [cell_id[1]+i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif in(cell_id[1], 2:M-1) && cell_id[3] == M
        for i in -1:1, j in [M-1, M, 1], k in [M-1, M, 1]
            key = [cell_id[1]+i, j, k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == 1 && in(cell_id[3], 2:M-1)
        for i in [M, 1, 2], j in [M-1, M, 1], k in -1:1
            key = [i, j, cell_id[3]+k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == M && in(cell_id[3], 2:M-1)
        for i in [M-1, M, 1], j in [M-1, M, 1], k in -1:1
            key = [i, j, cell_id[3]+k]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    end

    sum(k_t) > 0. ? v_t = sum(v_t) ./ sum(k_t) : v_t = zeros(3)

    return v_t
end

### ============== ### ============== ###
function get_neighbors(p_id, flock, box)
    if flock.p_cell_id[p_id] in box.bulk_cells
        # println("bulk_cells")
        flock.v_r[p_id] = check_bulk_cells(p_id, flock.p_cell_id[p_id], box.p_per_cell)
    elseif flock.p_cell_id[p_id] in box.bottom_cells
        # println("bottom_cells")
        flock.v_r[p_id] = check_bottom_cells(p_id, flock.p_cell_id[p_id], box.p_per_cell, box.L)
    elseif flock.p_cell_id[p_id] in box.top_cells
        # println("top_cells")
        flock.v_r[p_id] = check_top_cells(p_id, flock.p_cell_id[p_id], box.p_per_cell, box.L)
    elseif flock.p_cell_id[p_id] in box.left_wall
        # println("left_wall")
        flock.v_r[p_id] = check_left_cells(p_id, flock.p_cell_id[p_id], box.p_per_cell, box.L)
    elseif flock.p_cell_id[p_id] in box.right_wall
        # println("right_wall")
        flock.v_r[p_id] = check_right_cells(p_id, flock.p_cell_id[p_id], box.p_per_cell, box.L)
    elseif flock.p_cell_id[p_id] in box.front_wall
        # println("front_wall")
        flock.v_r[p_id] = check_front_cells(p_id, flock.p_cell_id[p_id], box.p_per_cell, box.L)
    elseif flock.p_cell_id[p_id] in box.back_wall
        # println("back_wall")
        flock.v_r[p_id] = check_back_cells(p_id, flock.p_cell_id[p_id], box.p_per_cell, box.L)
    elseif flock.p_cell_id[p_id] in box.corners
        # println("corners")
        flock.v_r[p_id] = check_corners(p_id, flock.p_cell_id[p_id], box.p_per_cell, box.L)
    end
end

### ============== ### ============== ###
function part_bulk_cell_interaction(p_id, flock, cell_parts)

    k = 0.
    v_r = zeros(Float64, 3)

    for j in cell_parts
        if norm(flock.pos[p_id] - flock.pos[j]) > 0. && norm(flock.pos[p_id] - flock.pos[j]) <= flock.r0
            k += 1.
            @. v_r = flock.vel[j]
        end
    end

    return k, v_r
end

### ============== ### ============== ###
function part_bound_cell_interaction(p_id, flock, cell_parts, L)

    k = 0.
    v_r = zeros(Float64, 3)
    d_v = zeros(Float64, 3)

    for j in cell_parts

        # compute relative vector  taking into account PBC
        for i in eachindex(flock.pos[p_id])
            δ = flock.pos[p_id][i] - flock.pos[j][i]
            δ > 0.5 * L ? d_v[i] = δ - L : d_v[i] = δ
        end

        # check if particles are within interaction range
        if  norm(d_v) > 0. && norm(d_v) <= flock.r0
            k += 1.
            @. v_r = flock.vel[j]
        end
    end

    return k, v_r
end

### ============== ### ============== ###
function update_part(pos, vel, v_r, η, L)

    q_r = Quaternion(zeros(Float64, 3))

    if norm(v_r) > zero(Float64)

        try
            # signal_angle = acos( dot(vel, v_r) / (norm(v_r) * norm(vel)) )
            signal_angle = acos( dot(vel, v_r) / norm(v_r) )

            # signal_angle = acos( dot(vel, normalize(v_r)) )
            q_r = qrotation(cross(vel, v_r), signal_angle + η * (2.0 * rand() * pi - pi)) * Quaternion(vel)
        catch error
            println("v_r")
            println(vel,"\t", v_r, "\t", dot(vel, v_r) / (norm(v_r) * norm(vel)))
        end

        # signal_angle = acos( dot(vel, v_r) / (norm(v_r) * norm(vel)) )
        # q_r = qrotation(cross(vel, v_r), signal_angle + η * (2.0 * rand() * pi - pi)) * Quaternion(vel)
    else

        noise = normalize!(randn(3))

        try
            signal_angle = acos( dot( noise, vel ) / norm(vel) )
            q_r = qrotation(cross(vel, noise), η * signal_angle ) * Quaternion(vel)
        catch error
            println("noise ")
            println(vel, "\t", noise, "\t", dot(noise, vel) / norm(vel))
        end

        # signal_angle = acos( dot( noise, vel ) / (norm(noise) * norm(vel)) )
        # q_r = qrotation(cross(vel, noise), η * signal_angle ) * Quaternion(vel)
    end

    ### ============== ###

    u_vel = normalize([q_r.v1, q_r.v2, q_r.v3])

    for i in eachindex(vel)

        vel[i] = u_vel[i]
        pos[i] += u_vel[i]

        if(pos[i] > L * 0.5) pos[i] -= L end
        if(pos[i] <= -L * 0.5) pos[i] += L end
    end

end

### ============== ### ============== ###

# u = [-0.650538, 0.702892, 0.287649]
# v = [-0.406586, 0.439308, 0.179781]
#
# norm(u)
# norm(v)
#
# acos(dot(u, v) / (norm(u) * norm(v)))
