### ============== ### ============== ###
### UTILITY FUNCTIONS FOR SVM 2D
### IN PERIODIC BOUNDARY CONDITION
### ============== ### ============== ###

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

    bulk_cells   :: Array{Array{Int64,1},2}
    bottom_cells :: Array{Array{Int64,1},1}
    top_cells    :: Array{Array{Int64,1},1}
    left_cells    :: Array{Array{Int64,1},1}
    right_cells   :: Array{Array{Int64,1},1}
    corners   :: Array{Array{Int64,1},2}

    p_per_cell :: Dict{Array{Int,1}, Array{Int, 1}}

    vel_cell :: Array{Array{Float64,1},2}
    center_cell :: Array{Array{Float64,1},2}

    function Box(L, M)

        ### ============== ### ============== ###
        ## indices of cells

        bulk_cells = [[i,j] for i in 2:M-1, j in 2:M-1]

        bottom_cells = [[i,M] for i in 2:M-1]
        top_cells = [[i,1] for i in 2:M-1]

        left_cells = [[1,j] for j in 2:M-1]
        right_cells = [[M,j] for j in 2:M-1]

        corners = [[i,j] for i in [1,M], j in [1,M]]

        ### ============== ### ============== ###

        # particles id's in each cell
        p_per_cell = Dict{Array{Int,1}, Array{Int, 1}}()

        # average velocity per cell
        vel_cell = Array{Array{Float64,1}}(M,M)

        # Initialize velocity vector per cell and compute center of each cell
        for i in 1:length(vel_cell)
            vel_cell[i] = zeros(Float64, 2)
        end

        # position of center of each box
        center_cell = Array{Array{Float64,1}}(M,M)

        r = linspace(0., L, M+1)
        mp = collect(r[1:length(r) - 1] + 0.5 * step(r))

        for i in 1:M, j in 1:M
            # println(i,"\t", j, "\t", k)
            center_cell[i,j] = [mp[i], mp[j]]
        end

        new(L, float(M), bulk_cells, bottom_cells, top_cells, left_cells, right_cells, corners, p_per_cell, vel_cell, center_cell)
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

* p_cell_id -> particle's cell id
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
        pos = [ [L*rand(), L*rand()] for i in 1:N ]
        vel = v0 * [ normalize([2*rand() - 1, 2*rand() - 1]) for i in 1:N ]

        # interaction vector
        v_r = [zeros(2) for i in 1:N]

        # connectivity per particle
        k_r = zeros(Float64, N)

        # box index for each particle in each dimension
        p_cell_id = [zeros(Int,2) for i in 1:N]

        new(N, v0, r0, dt, η, pos, vel, v_r, k_r, p_cell_id)
    end
end

## ============== ### ============== ###
function av_vel_cell(box, vel)
    for k in keys(box.p_per_cell)
        box.vel_cell[k[1], k[2]] = mean([vel[i] for i in box.p_per_cell[k]])
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

    k_t = Vector{Float64}(9)
    v_t = Vector{Vector{Float64}}(9)

    c = 1

    # println("center cell: ", cell_id)

    for i in -1:1, j in -1:1
        key = [cell_id[1]+i, cell_id[2]+j]
        k_t[c], v_t[c] = part_bulk_cell_interaction(p_id, flock, get!(p_per_cell, key, []))
        c += 1
        # println(key, "\t", get!(p_per_cell, key, []))
    end

    sum(k_t) > 0. ? v_t = sum(v_t) ./ sum(k_t) : v_t = zeros(2)

    return v_t
end

### ============== ### ============== ###
function check_corners(p_id, cell_id, p_per_cell, L)

    k_t = Vector{Float64}(9)
    v_t = Vector{Vector{Float64}}(9)

    c = 1
    # println("corner cell: ", cell_id)

    if cell_id[1] == 1 && cell_id[2] == 1
        for i in [M, 1, 2], j in [M, 1, 2]
            key = [i, j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == 1 && cell_id[2] == M
        for i in [M, 1, 2], j in [M-1, M, 1]
            key = [i, j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == M && cell_id[2] == 1
        for i in [M-1, M, 1], j in [M, 1, 2]
            key = [i, j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == M && cell_id[2] == M
        for i in [M-1, M, 1], j in [M-1, M, 1]
            key = [i, j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    end

    sum(k_t) > 0. ? v_t = sum(v_t) ./ sum(k_t) : v_t = zeros(2)

    return v_t
end

### ============== ### ============== ###
function check_bottom_cells(p_id, cell_id, p_per_cell, L)

    k_t = Vector{Float64}(9)
    v_t = Vector{Vector{Float64}}(9)

    c = 1
    # println("center cell: ", cell_id)

    # if in(cell_id[1], 2:M-1) && cell_id[2] == M
        for i in -1:1, j in [M-1, M, 1]
            key = [cell_id[1] + i, j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    # end

    sum(k_t) > 0. ? v_t = sum(v_t) ./ sum(k_t) : v_t = zeros(2)

    return v_t
end

### ============== ### ============== ###
function check_top_cells(p_id, cell_id, p_per_cell, L)
    k_t = Vector{Float64}(9)
    v_t = Vector{Vector{Float64}}(9)

    c = 1
    # println("center cell: ", cell_id)

    # if in(cell_id[1], 2:M-1) && cell_id[2] == 1
        for i in -1:1, j in [M, 1, 2]
            key = [cell_id[1] + i, j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    # end

    sum(k_t) > 0. ? v_t = sum(v_t) ./ sum(k_t) : v_t = zeros(2)

    return v_t
end

### ============== ### ============== ###
function check_left_cells(p_id, cell_id, p_per_cell, L)

    k_t = Vector{Float64}(9)
    v_t = Vector{Vector{Float64}}(9)

    c = 1
    # println("center cell: ", cell_id)

    # if cell_id[1] == 1 && in(cell_id[2], 2:M-1)
        for i in [M, 1, 2], j in -1:1
            key = [i, cell_id[2] + j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    # end

    sum(k_t) > 0. ? v_t = sum(v_t) ./ sum(k_t) : v_t = zeros(2)

    return v_t
end

### ============== ### ============== ###
function check_right_cells(p_id, cell_id, p_per_cell, L)

    k_t = Vector{Float64}(9)
    v_t = Vector{Vector{Float64}}(9)

    c = 1
    # println("center cell: ", cell_id)

    # if cell_id[1] == M && in(cell_id[2], 2:M-1)
        for i in [M-1, M, 1], j in -1:1
            key = [i, cell_id[2] + j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    # end

    sum(k_t) > 0. ? v_t = sum(v_t) ./ sum(k_t) : v_t = zeros(2)

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
    elseif flock.p_cell_id[p_id] in box.left_cells
        # println("left_cells")
        flock.v_r[p_id] = check_left_cells(p_id, flock.p_cell_id[p_id], box.p_per_cell, box.L)
    elseif flock.p_cell_id[p_id] in box.right_cells
        # println("right_cells")
        flock.v_r[p_id] = check_right_cells(p_id, flock.p_cell_id[p_id], box.p_per_cell, box.L)
    elseif flock.p_cell_id[p_id] in box.corners
        # println("corners")
        flock.v_r[p_id] = check_corners(p_id, flock.p_cell_id[p_id], box.p_per_cell, box.L)
    end
end

### ============== ### ============== ###
function part_bulk_cell_interaction(p_id, flock, cell_parts)

    k = 0.
    v_r = zeros(Float64, 2)

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
    v_r = zeros(Float64, 2)
    d_v = zeros(Float64, 2)

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

    prop_angle = atan2(vel[2], vel[1])
    signal_angle = 0.0
    u_vel = zeros(Float64, 2)

    v_r[1] != 0.0 || v_r[2] != 0.0 ? signal_angle = atan2(v_r[2], v_r[1]) - prop_angle : signal_angle = 0.0

    total_angle = signal_angle + η * (2.0 * rand() * pi - pi);

    c = cos(total_angle)
    s = sin(total_angle)

    u_vel[1] = vel[1]*c - vel[2]*s;
    u_vel[2] = vel[1]*s + vel[2]*c;

    for i in eachindex(vel)

        vel[i] = u_vel[i]
        pos[i] += u_vel[i]

        # if(pos[i] > L * 0.5) pos[i] -= L end
        # if(pos[i] <= -L * 0.5) pos[i] += L end

        if(pos[i] > L ) pos[i] -= L end
        if(pos[i] <= 0.0) pos[i] += L end
    end

end

### ============== ### ============== ###
