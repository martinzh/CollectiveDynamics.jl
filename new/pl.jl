####========================================####
####           TYPE DEFINITIONS             ####
####========================================####

type Param
    M::Int64
    N::Int64
    v0::Float64
    r0::Float64
    eta::Float64
    omega::Float64
    bound::Float64

    # Inicializa todo a ceros
    Param() = new(0, 0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

####========================================####

type Flock

    pos::SharedArray
    vel::SharedArray
    v_r::SharedArray
    v_n::SharedArray
    rij::SharedArray
    noise::SharedArray
    poski::SharedArray
    nij::SharedArray
    rij_ids::SharedArray

    #Inicializa
    Flock(N::Int64, M::Int64) = new(
        convert( SharedArray, zeros(2N) ), # pos
        convert( SharedArray, zeros(2N) ), # vel
        convert( SharedArray, zeros(2N) ), # v_r
        convert( SharedArray, zeros(2N) ), # v_n
        convert( SharedArray, zeros(N, N) ),# rij
        convert( SharedArray, zeros(N) ),  # noise
        convert( SharedArray, zeros(N + M) ), # poski
        convert( SharedArray, zeros(M) ), # nij
        convert(SharedArray, Array{Int64}(div(N * (N-1), 2), 2)) #rij id's
    )

end

####========================================####
####      SHARED ARRAYS UILITY FUNCTIONS    ####
####========================================####

function init_vels(N::Int64, v0::Float64, vels::SharedArray)

    for i in 1:2:2N
        vx = -1.0 + 2.0 * rand()
        vy = -1.0 + 2.0 * rand()

        # v = norm(Float64[vx,vy])
        v = sqrt(vx^2 + vy^2)

        vels[i]     = v0 * (vx / v)
        vels[i + 1] = v0 * (vy / v)
    end

end

function init_system(flock::Flock, param::Param)

    # posiciones iniciales
    # particulas se inicializan en una caja de tamanio w
    flock.pos = Float64[ -param.bound + rand() * 2.0 * param.bound for i in 1:(2 * param.N) ]

    # velocidades
    # velocidad inicial con direccion aleatoria y norma v0
    init_vels(param.N, param.v0, flock.vel)

    #inicializa red de interaccion no local
    nij, poski  = make_IN(param.N, param.M)
    flock.nij   = convert(SharedArray, nij)
    flock.poski = convert(SharedArray, poski)

    assign_ids(flock.rij_ids, param.N)

end
###========================================###

function make_dir(path::ASCIIString)
    try
        run(`mkdir $path`)
    catch y
        println(typeof(y))
    end
end

####========================================####

# Range assignation for 1D sh array
@everywhere function myrange(q::SharedArray)
    idx = indexpids(q)
    if idx == 0
        # This worker is not assigned a piece
        return 1:0, 1:0
    end
    nchunks = length(procs(q))
    splits = [round(Int, s) for s in linspace(0, size(q,1), nchunks+1)]
    # println(splits[idx]+1:splits[idx+1])
    # println(myid())
    return splits[idx]+1:splits[idx+1]
end

####===============================####

@everywhere function my_rij_range(N::Int64)
    idx = myid()-1
    if idx == 0
        # This worker is not assigned a piece
        return 1:0, 1:0
    end
    nchunks = length(procs()[2:end])
    splits = [round(Int, s) for s in linspace(0, N, nchunks+1)]
    # println(splits[idx]+1:splits[idx+1], ", ", typeof(splits[idx]+1:splits[idx+1]) )
    return splits[idx]+1:splits[idx+1]
end

####===============================####

function calc_ranges(N::Int64)
    nchunks = length(procs()[2:end])
    splits = [round(Int, s) for s in linspace(0, N, nchunks+1)]
    return UnitRange[splits[idx-1]+1:splits[idx] for idx in procs()[2:end]]
end
####===============================####

function assign_ids(pair_dists::SharedArray, N::Int64)

    k = 1

    for i in 1:N, j in (i+1):N
        pair_dists[k,1] = i
        pair_dists[k,2] = j
        k += 1
    end

end

####===============================####

# @everywhere function rot_move_chunk(pos::SharedArray, vel::SharedArray, v_r::SharedArray, v_n::SharedArray, noise::SharedArray, dt::Float64, omega::Float64, eta::Float64, id_range::UnitRange)
@everywhere function rot_move_chunk(pos::SharedArray, vel::SharedArray, v_r::SharedArray, v_n::SharedArray, dt::Float64, omega::Float64, eta::Float64, id_range::UnitRange, vx::Float64, vy::Float64)

    @inbounds for i in id_range

        # prop_angle = atan2(vel[2i], vel[2i - 1])

        # i_vx = v_r[2i - 1]
        # i_vy = v_r[2i]

        # if i_vx != 0.0 || i_vy != 0.0
        #     loc_angle = atan2(i_vy, i_vx) - prop_angle
        # else
        #     loc_angle = 0.0
        # end

        # i_vx != 0.0 || i_vy != 0.0 ? loc_angle = atan2(i_vy, i_vx) - prop_angle : loc_angle = 0.0

        # tot_angle = omega * loc_angle + (noise[i] * 2.0 * pi - pi) * eta

        # i_vx = v_n[2i - 1]
        # i_vy = v_n[2i]

        # if i_vx != 0.0 || i_vy != 0.0
        #     nonloc_angle = atan2(i_vy, i_vx) - prop_angle
        # else
        #     nonloc_angle = 0.0
        # end

        # i_vx != 0.0 || i_vy != 0.0 ? nonloc_angle = atan2(i_vy, i_vx) - prop_angle : nonloc_angle = 0.0

        # loc_angle    = atan2(v_r[2i - 1], v_r[2i]) - prop_angle
        # nonloc_angle = atan2(v_n[2i - 1], v_n[2i]) - prop_angle

        # tot_angle = omega * loc_angle + (1.0 - omega) * nonloc_angle + (noise[i] * 2.0 * pi - pi) * eta

        # tot_angle = omega * loc_angle + (1.0 - omega) * nonloc_angle + (rand() * 2.0 * pi - pi) * eta

        tot_angle = omega * (atan2(v_r[2i - 1], v_r[2i]) - atan2(vel[2i], vel[2i - 1])) + (1.0 - omega) * (atan2(v_n[2i - 1], v_n[2i]) - atan2(vel[2i - 1], vel[2i])) + (rand() * 2.0 * pi - pi) * eta

        vx = vel[2i - 1]*cos(tot_angle) - vel[2i]*sin(tot_angle)
        vy = vel[2i - 1]*sin(tot_angle) + vel[2i]*cos(tot_angle)

        vel[2i - 1] = vx
        vel[2i]     = vy

        pos[2i - 1] += vel[2i - 1] * dt
        pos[2i]     += vel[2i] * dt

    end

end

@everywhere rot_move(pos::SharedArray, vel::SharedArray, v_r::SharedArray, v_n::SharedArray, noise::SharedArray, dt::Float64, omega::Float64, eta::Float64, N::Int64) = rot_move_chunk(pos, vel, v_r, v_n, noise, dt, omega, eta, my_rij_range(N))

####===============================####

@everywhere function rot_move_chunk_bound(pos::SharedArray, vel::SharedArray, v_r::SharedArray, v_n::SharedArray, noise::SharedArray, dt::Float64, omega::Float64, eta::Float64, bound::Float64, id_range::UnitRange)

    @inbounds for i in id_range

        prop_angle = atan2(vel[2i], vel[2i - 1])

        i_vx = v_r[2i - 1]
        i_vy = v_r[2i]

        if i_vx != 0.0 || i_vy != 0.0
            loc_angle = atan2(i_vy, i_vx) - prop_angle
        else
            loc_angle = 0.0
        end

        # tot_angle = omega * loc_angle + (noise[i] * 2.0 * pi - pi) * eta

        i_vx = v_n[2i - 1]
        i_vy = v_n[2i]

        if i_vx != 0.0 || i_vy != 0.0
            nonloc_angle = atan2(i_vy, i_vx) - prop_angle
        else
            nonloc_angle = 0.0
        end

        tot_angle = omega * loc_angle + (1.0 - omega) * nonloc_angle + (noise[i] * 2.0 * pi - pi) * eta
        tot_angle = omega * loc_angle + (1.0 - omega) * nonloc_angle + (noise[i] * 2.0 * pi - pi) * eta

        vx = vel[2i - 1]*cos(tot_angle) - vel[2i]*sin(tot_angle)
        vy = vel[2i - 1]*sin(tot_angle) + vel[2i]*cos(tot_angle)

        vel[2i - 1] = vx
        vel[2i]     = vy

        pos[2i - 1] += vel[2i - 1] * dt
        pos[2i]     += vel[2i] * dt

        #  periodic bounds
        if pos[2i - 1] > (bound * 0.5)
            pos[2i - 1] -= bound
        elseif pos[2i - 1] <= (-bound * 0.5)
            pos[2i - 1] += bound
        end

        if pos[2i] > (bound * 0.5)
            pos[2i] -= bound
        elseif pos[2i] <= (-bound * 0.5)
            pos[2i] += bound
        end

    end

end

@everywhere rot_move_bound(pos::SharedArray, vel::SharedArray, v_r::SharedArray, v_n::SharedArray, noise::SharedArray, dt::Float64, omega::Float64, eta::Float64, bound::Float64, N::Int64) = rot_move_chunk_bound(pos, vel, v_r, v_n, noise, dt, omega, eta, bound, my_rij_range(N))

####===============================####

@everywhere function calc_Rij_chunk(pos::SharedArray, rij::SharedArray, N::Int64, r0::Float64, id_range::UnitRange)

    # r02 = r0*r0

    @inbounds for i in id_range, j in (i+1):N

        # d =  (pos[2i - 1] - pos[2j - 1])^2 + (pos[2i] - pos[2j])^2

        # if d <= r0*r0
        #     # rij[j,i] = 1.0
        #     rij[j,i] = rij[i,j] = 1.0
        # else
        #     # rij[j,i] = 0.0
        #     rij[j,i] = rij[i,j] = 0.0
        # end

        # rij[j, i] = rij[i, j] = 666.0

        (pos[2i - 1] - pos[2j - 1])^2 + (pos[2i] - pos[2j])^2 <= r0*r0 ? rij[j,i] = rij[i,j] = 1.0 : rij[j,i] = rij[i,j] = 0.0

    end

end

# WRAPPER
@everywhere calc_Rij(pos::SharedArray, rij::SharedArray, N::Int64, r0::Float64) = calc_Rij_chunk(pos, rij, N, r0, my_rij_range(N))

####===============================####

@everywhere function calc_Rij_chunk_id(pos::SharedArray, rij::SharedArray, rij_ids::SharedArray, r0::Float64, id_range::UnitRange)

    @inbounds @simd for k in id_range

        # i = rij_ids[k, 1]
        # j = rij_ids[k, 2]

        # d =  (pos[2i - 1] - pos[2j - 1])^2 + (pos[2i] - pos[2j])^2
        #
        # if d <= r0 * r0
        #     # rij[j,i] = 1.0
        #     rij[j,i] = rij[i,j] = 1.0
        # else
        #     # rij[j,i] = 0.0
        #     rij[j,i] = rij[i,j] = 0.0
        # end

        # (pos[2i - 1] - pos[2j - 1])^2 + (pos[2i] - pos[2j])^2 <= r0*r0 ? rij[j,i] = rij[i,j] = 1.0 : rij[j,i] = rij[i,j] = 0.0

        (pos[2*rij_ids[k, 1] - 1] - pos[2*rij_ids[k, 2] - 1])^2 + (pos[2*rij_ids[k, 1]] - pos[2*rij_ids[k, 2]])^2 <= r0*r0 ? rij[rij_ids[k, 2], rij_ids[k, 1]] = rij[rij_ids[k, 1], rij_ids[k, 2]] = 1.0 : rij[rij_ids[k, 2], rij_ids[k, 1]] = rij[rij_ids[k, 1], rij_ids[k, 2]] = 0.0

        # (pos[2*rij_ids[k, 1] - 1] - pos[2*rij_ids[k, 2] - 1])^2 + (pos[2*rij_ids[k, 1]] - pos[2*rij_ids[k, 2]])^2 <= r0*r0 ? rij[rij_ids[k, 2], rij_ids[k, 1]] = 1.0 : rij[rij_ids[k, 2], rij_ids[k, 1]] = 0.0

        # rij[j,i] = 666.0
    end

end

# WRAPPER
@everywhere calc_Rij_id(pos::SharedArray, rij::SharedArray, rij_ids::Int64, r0::Float64, tot_dists::Int64) = calc_Rij_chunk(pos, rij, rij_ids, r0, my_rij_range(tot_dists))

####===============================####

function dists(flock::Flock, param::Param)
    @sync begin
        for p in procs()[2:end]
            @async remotecall_wait(p, calc_Rij, flock.pos, flock.rij, param.N, param.r0)
        end
    end
end

function dists_id(flock::Flock, param::Param, rij_ranges::Array{UnitRange})
    @sync begin
        for p in procs()[2:end]
            @async remotecall_wait(p, calc_Rij_chunk_id, flock.pos, flock.rij, flock.rij_ids, param.r0, rij_ranges[p - 1])
        end
    end
end

####===============================####

@everywhere function loc_vels_chunk(vel::SharedArray, v_r::SharedArray, rij::SharedArray, N::Int64, id_range::UnitRange)
# @everywhere function loc_vels_chunk(vel::SharedArray, v_r::SharedArray, rij::SharedArray, N::Int64, id_range::UnitRange, vx::Float64, vy::Float64)

    @inbounds for i in id_range

        vx = vel[2i - 1]
        vy = vel[2i]

        ki = 1.0 + sum(rij[: ,i])

        @inbounds for j in 1:N

            # adj = rij[ j, i ]
            # vx += vel[2j - 1] * adj
            # vy += vel[2j]     * adj

            vx += vel[2j - 1] * rij[ j, i ]
            vy += vel[2j]     * rij[ j, i ]
        end

        # for j in findn(rij[:, i])
        #
        #     vx += vel[2j - 1]
        #     vy += vel[2j]
        #
        # end

        v_r[2i - 1] = vx / ki
        v_r[2i]     = vy / ki
    end

end

@everywhere loc_vels(vel::SharedArray, v_r::SharedArray, rij::SharedArray, N::Int64) = loc_vels_chunk(vel, v_r, rij, N, my_rij_range(N))

####===============================####
@everywhere function non_loc_vels_chunk(vel::SharedArray, v_n::SharedArray, nij::SharedArray, poski::SharedArray, id_range::UnitRange)

    @inbounds for i in id_range

        vx = 0.0
        vy = 0.0

        ki = nij[ poski[i] ]

        if ki > 0

            @inbounds for j in 1:ki

                # k = nij[ poski[i] + j ]
                #
                # vx += vel[2k - 1]
                # vy += vel[2k]

                vx += vel[2 * nij[ poski[i] + j ] - 1]
                vy += vel[2 * nij[ poski[i] + j ]]

            end

            v_n[2i - 1] = vx / convert(Float64, ki)
            v_n[2i]     = vy / convert(Float64, ki)

        else
            v_n[2i - 1] = 0.0
            v_n[2i]     = 0.0
        end

    end

end

@everywhere non_loc_vels(vel::SharedArray, v_n::SharedArray, nij::SharedArray, poski::SharedArray, N::Int64) = non_loc_vels_chunk(vel, v_n, nij, poski, my_rij_range(N))

####===============================####

function calc_vels(flock::Flock, param::Param)

    @sync begin
        for p in procs(pos)
            @async remotecall_wait(p, loc_vels, flock.vel, flock.v_r, flock.rij, param.N)
            @async remotecall_wait(p, non_loc_vels, flock.vel, flock.v_n, flock.nij, flock.poski, param.N)
        end
    end

end

####===============================####

function evol(flock::Flock, param::Param, vals::Array{Float64,2}, T::Int64, dt::Float64)
    for t in 1:T

        @sync begin
            for p in procs(flock.pos)
                @async remotecall_wait(p, calc_Rij, flock.pos, flock.rij, param.N, param.r0)
            end
        end

        @sync begin
            for p in procs(flock.pos)
                @async remotecall_wait(p, loc_vels, flock.vel, flock.v_r, flock.rij, param.N)
                @async remotecall_wait(p, non_loc_vels, flock.vel, flock.v_n, flock.nij, flock.poski, param.N)
            end
        end

        # flock.noise = convert(SharedArray, rand(param.N))
        flock.noise = rand(param.N)

        @sync begin
            for p in procs(flock.pos)
                @async remotecall_wait(p, rot_move, flock.pos, flock.vel, flock.v_r, flock.v_n, flock.noise, dt, param.omega, param.eta, param.N)
            end
        end

        vals[:,t] = flock.pos.s
    end
end

####===============================####

function evol_bound(flock::Flock, param::Param, vals::Array{Float64,2}, T::Int64, dt::Float64)
    for t in 1:T

        @sync begin
            for p in procs(flock.pos)
                @async remotecall_wait(p, calc_Rij, flock.pos, flock.rij, param.N, param.r0)
            end
        end

        @sync begin
            for p in procs(flock.pos)
                @async remotecall_wait(p, loc_vels, flock.vel, flock.v_r, flock.rij, param.N)
                @async remotecall_wait(p, non_loc_vels, flock.vel, flock.v_n, flock.nij, flock.poski, param.N)
            end
        end

        flock.noise = rand(param.N)

        @sync begin
            for p in procs(flock.pos)
                @async remotecall_wait(p, rot_move_bound, flock.pos, flock.vel, flock.v_r, flock.v_n, flock.noise, dt, param.omega, param.eta, param.bound, param.N)
            end
        end

        vals[:,t] = flock.pos.s
    end
end

####===============================####

function evol_range(flock::Flock, param::Param, vals::Array{Float64,2}, T::Int64, dt::Float64, range::Array{UnitRange}, rij_range::Array{UnitRange})
    for t in 1:T

        @sync begin
            for p in procs(flock.pos)
                # @async remotecall_wait(p, calc_Rij_chunk, flock.pos, flock.rij, param.N, param.r0, range[p-1])
                @async remotecall_wait(p, calc_Rij_chunk_id, flock.pos, flock.rij, flock.rij_ids, param.r0, rij_range[p-1])
            end
        end

        @sync begin
            for p in procs(flock.pos)
                @async remotecall_wait(p, loc_vels_chunk, flock.vel, flock.v_r, flock.rij, param.N, range[p-1])
                @async remotecall_wait(p, non_loc_vels_chunk, flock.vel, flock.v_n, flock.nij, flock.poski, range[p-1])
            end
        end

        flock.noise = rand(param.N)

        @sync begin
            for p in procs(flock.pos)
                @async remotecall_wait(p, rot_move_chunk, flock.pos, flock.vel, flock.v_r, flock.v_n, flock.noise, dt, param.omega, param.eta, range[p-1])
            end
        end

        vals[:,t] = flock.pos.s
    end
end

####===============================####

function evol_step_range(flock::Flock, param::Param, dt::Float64, range::Array{UnitRange}, rij_range::Array{UnitRange}, vx::Array{Float64,1}, vy::Array{Float64,1})

    @sync begin
        for p in procs(flock.pos)
            @async remotecall_wait(p, calc_Rij_chunk_id, flock.pos, flock.rij, flock.rij_ids, param.r0, rij_range[p-1])
        end
    end

    @sync begin
        for p in procs(flock.pos)
            @async remotecall_wait(p, loc_vels_chunk, flock.vel, flock.v_r, flock.rij, param.N, range[p-1])
            # @async remotecall_wait(p, loc_vels_chunk, flock.vel, flock.v_r, flock.rij, param.N, range[p-1], vx[p-1], vy[p-1])
            @async remotecall_wait(p, non_loc_vels_chunk, flock.vel, flock.v_n, flock.nij, flock.poski, range[p-1])
        end
    end

    # flock.noise = rand(param.N)

    @sync begin
        for p in procs(flock.pos)
            # @async remotecall_wait(p, rot_move_chunk, flock.pos, flock.vel, flock.v_r, flock.v_n, flock.noise, dt, param.omega, param.eta, range[p-1])
            @async remotecall_wait(p, rot_move_chunk, flock.pos, flock.vel, flock.v_r, flock.v_n, dt, param.omega, param.eta, range[p-1], vx[p-1], vy[p-1])
        end
    end

end

####===============================####

function evol_bound_range(flock::Flock, param::Param, vals::Array{Float64,2}, T::Int64, dt::Float64, range::Array{UnitRange})
    for t in 1:T

        @sync begin
            for p in procs(flock.pos)
                @async remotecall_wait(p, calc_Rij_chunk, flock.pos, flock.rij, param.r0, range[p-1])
            end
        end

        @sync begin
            for p in procs(flock.pos)
                @async remotecall_wait(p, loc_vels_chunk, flock.vel, flock.v_r, flock.rij, range[p-1])
                @async remotecall_wait(p, non_loc_vels_chunk, flock.vel, flock.v_n, flock.nij, flock.poski, range[p-1])
            end
        end

        flock.noise = rand(param.N)

        @sync begin
            for p in procs(flock.pos)
                @async remotecall_wait(p, rot_move_bound_chunk, flock.pos, flock.vel, flock.v_r, flock.v_n, flock.noise, dt, param.omega, param.eta, param.bound, range[p-1])
            end
        end

        vals[:,t] = flock.pos.s
    end
end

####===============================####
####        TEST FUNCTIONS         ####
####===============================####

function print_range(q)
    @sync begin
        for p in procs(q)
            @async remotecall_wait(p, myrange, q)
        end
    end
end

function print_rij_range(N)
    @sync begin
        for p in procs()
            @async remotecall_wait(p, my_rij_range, N)
        end
    end
end

###========================================###
### FUNCTIONS FOR CREATING INTERACTION NET ###
###========================================###

###========================================###

function checkId(id::Int64, n::Int64)
    for i in 0:( n-1 )
        if in(id, ( i * n+1 ):( n * (i+1) ) )
            return i+1
        end
    end
end

###========================================####

function setLinks!(m::Int64, n::Int64, numLinks::Array{Int64, 1}, adj_vec::Array{Int64, 1})
    for i in 1:m
        while true
            link = rand(1:n*n)
            # link = int(rand(rngLinks)*n*n)+1
            if adj_vec[link] < 1
                adj_vec[link] = 1
                numLinks[checkId(link, n)] += 1
                break
            end
        end
    end
end

###========================================####

function reArrange(m::Int64, n::Int64, numLinks::Array{Int64, 1}, adj_vec::Array{Int64, 1})

    Nij = zeros(Int64, m+n)
    posNumLinks = zeros(Int64, n)

    k = 1
    for i in 1:n
        Nij[k] = numLinks[i]
        posNumLinks[i] = k
        k += numLinks[i] + 1
    end

    links = find(x -> x==1, adj_vec)

    k = 0
    for i in 1:n
        for j in 1:numLinks[i]
            part = links[j+k] % n
            if part != 0 # part - 1 -> [0,n-1] for GPU
                Nij[posNumLinks[i]+j] = part
            else
                Nij[posNumLinks[i]+j] = n
            end
        end
        k += numLinks[i]
    end

    return Nij, posNumLinks
end

###========================================####

function make_IN(n::Int64, m::Int64)
    adj_vec = reshape(2*eye(Int64, n), n*n)
    numLinks = zeros(Int64, n)

    setLinks!(m, n, numLinks, adj_vec)
    return reArrange(m, n, numLinks, adj_vec)
end
