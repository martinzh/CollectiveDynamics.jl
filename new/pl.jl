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

    #Inicializa
    Flock(N::Int64, M::Int64) = new(
        convert(SharedArray, zeros(2N)), # pos
        convert(SharedArray, zeros(2N)), # vel
        convert(SharedArray, zeros(2N)), # v_r
        convert(SharedArray, zeros(2N)), # v_n
        convert(SharedArray, zeros(N,N)),# rij
        convert(SharedArray, zeros(N)),  # noise
        convert(SharedArray, zeros(N + M)), # poski
        convert(SharedArray, zeros(M)) # nij
    )

end

####========================================####
####      SHARED ARRAYS UILITY FUNCTIONS    ####
####========================================####

function init_vels(N::Int64, v0::Float64, vels::SharedArray)

    for i in 1:2:2N
        vx = -1.0 + 2.0 * rand()
        vy = -1.0 + 2.0 * rand()

        v = norm(Float64[vx,vy])

        vels[i]     = v0 * (vx / v)
        vels[i + 1] = v0 * (vy / v)
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
    splits = [round(Int, s) for s in linspace(0,size(q,1),nchunks+1)]
    # println(splits[idx]+1:splits[idx+1])
    # println(myid())
    splits[idx]+1:splits[idx+1]
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
    splits[idx]+1:splits[idx+1]
end

####===============================####

@everywhere function rot_move_chunk(pos::SharedArray, vel::SharedArray, v_r::SharedArray, noise::SharedArray, dt::Float64, omega::Float64, eta::Float64, id_range::UnitRange)

    for i in id_range

        prop_angle = atan2(vel[2i], vel[2i - 1])

        i_vx = v_r[2i - 1]
        i_vy = v_r[2i]

        if i_vx != 0.0 || i_vy != 0.0
            loc_angle = atan2(i_vy, i_vx) - prop_angle
        else
            loc_angle = 0.0
        end

        # tot_angle = omega * loc_angle + (1.0 - omega) * nonloc_angle + (noise[i] * 2.0 * pi - pi) * eta

        tot_angle = omega * loc_angle + (noise[i] * 2.0 * pi - pi) * eta

        # i_vx = nl_vel[2i - 1]
        # i_vy = nl_vel[2i]
        #
        # if i_vx != 0.0 || i_vy != 0.0
        #     nonloc_angle = atan2(i_vy, i_vx) - prop_angle
        # else
        #     nonloc_angle = 0.0
        # end

        # tot_angle = omega * loc_angle + (1.0 - omega) * nonloc_angle + (noise[i] * 2.0 * pi - pi) * eta

        vx = vel[2i - 1]*cos(tot_angle) - vel[2i]*sin(tot_angle)
        vy = vel[2i - 1]*sin(tot_angle) + vel[2i]*cos(tot_angle)

        vel[2i - 1] = vx
        vel[2i]     = vy

        pos[2i - 1] += vel[2i - 1] * dt
        pos[2i]     += vel[2i] * dt

    end

end

@everywhere rot_move(pos::SharedArray, vel::SharedArray, v_r::SharedArray, noise::SharedArray, dt::Float64, omega::Float64, eta::Float64, N::Int64) = rot_move_chunk(pos, vel, v_r, noise, dt, omega::Float64, eta::Float64, my_rij_range(N))

####===============================####

@everywhere function rot_move_chunk_bound(pos::SharedArray, vel::SharedArray, v_r::SharedArray, noise::SharedArray, dt::Float64, omega::Float64, eta::Float64, bound::Float64, id_range::UnitRange)

    for i in id_range

        prop_angle = atan2(vel[2i], vel[2i - 1])

        i_vx = v_r[2i - 1]
        i_vy = v_r[2i]

        if i_vx != 0.0 || i_vy != 0.0
            loc_angle = atan2(i_vy, i_vx) - prop_angle
        else
            loc_angle = 0.0
        end

        tot_angle = omega * loc_angle + (noise[i] * 2.0 * pi - pi) * eta

        # i_vx = nl_vel[2i - 1]
        # i_vy = nl_vel[2i]
        #
        # if i_vx != 0.0 || i_vy != 0.0
        #     nonloc_angle = atan2(i_vy, i_vx) - prop_angle
        # else
        #     nonloc_angle = 0.0
        # end

        # tot_angle = omega * loc_angle + (1.0 - omega) * nonloc_angle + (noise[i] * 2.0 * pi - pi) * eta

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

@everywhere rot_move_bound(pos::SharedArray, vel::SharedArray, v_r::SharedArray, noise::SharedArray, dt::Float64, omega::Float64, eta::Float64, bound::Float64, N::Int64) = rot_move_chunk_bound(pos, vel, v_r, noise, dt, omega::Float64, eta::Float64, bound::Float64, my_rij_range(N))

####===============================####

@everywhere function calc_Rij_chunk(pos::SharedArray, rij::SharedArray, N::Int64, r0::Float64, id_range::UnitRange)

    for i in id_range, j in (i+1):N

        d =  (pos[2i - 1] - pos[2j - 1])^2 + (pos[2i] - pos[2j])^2

        if d <= r0 * r0
            # rij[j,i] = 1
            rij[j,i] = rij[i,j] = 1.0
        else
            rij[j,i] = rij[i,j] = 0.0
        end
    end

end

# WRAPPER
@everywhere calc_Rij(pos::SharedArray, rij::SharedArray, N::Int64, r0::Float64) = calc_Rij_chunk(pos, rij, N, r0, my_rij_range(N))

####===============================####

# function dists(pos::SharedArray, rij::SharedArray, N::Int64, r0::Float64)
function dists(flock::Flock, param::Param)
    @sync begin
        for p in procs()[2:end]
            # @async remotecall_wait(p, calc_Rij, pos, rij, N, r0)
            @async remotecall_wait(p, calc_Rij, flock.pos, flock.rij, param.N, param.r0)
        end
    end
end

####===============================####

@everywhere function loc_vels_chunk(vel::SharedArray, v_r::SharedArray, rij::SharedArray, N::Int64, id_range::UnitRange)

    for i in id_range

        vx = vel[2i - 1]
        vy = vel[2i]

        ki = 1.0 + sum(rij.s[: ,i])

        for j in 1:N
            vx += vel[2j - 1] * rij[j,i]
            vy += vel[2j]     * rij[j,i]
        end

        v_r[2i - 1] = vx / ki
        v_r[2i]     = vy / ki
    end

end

@everywhere loc_vels(vel::SharedArray, v_r::SharedArray, rij::SharedArray, N::Int64) = loc_vels_chunk(vel, v_r, rij, N, my_rij_range(N))

####===============================####

function calc_vels(flock::Flock, param::Param)

    @sync begin
        for p in procs(pos)
            @async remotecall_wait(p, loc_vels, flock.vel, flock.v_r, flock.rij, param.N)
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
            end
        end

        @sync begin
            for p in procs(flock.pos)
                @async remotecall_wait(p, rot_move, flock.pos, flock.vel, flock.v_r, flock.noise, dt, param.omega, param.eta, param.N)
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
            end
        end

        @sync begin
            for p in procs(flock.pos)
                @async remotecall_wait(p, rot_move_bound, flock.pos, flock.vel, flock.v_r, flock.noise, dt, param.omega, param.eta, param.bound, param.N)
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

####===============================####
