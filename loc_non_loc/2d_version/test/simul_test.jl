# Modelo movimiento colectivo
# Paralelo con SharedArrays

# using CollectiveDynamics
addprocs(4)

@everywhere include("/Users/martinC3/GitRepos/CollectiveDynamics.jl/new/pl.jl")

using ProfileView

using PyPlot

####========================================####
####             PARAMETERS                 ####
####========================================####

param = Param()

steps = frec = reps = t = 0
rho = m_frac = w = 0.0
dt = 1.0
l = 0.1

param.v0 = 1.0

param.r0 = param.v0 * dt / l
param.omega = 0.5 # peso interacciones

if length(ARGS) == 7

    #  Parametros de linea de comando

    param.N   = parse(Int64,   ARGS[1]);
    m_frac    = parse(Float64, ARGS[2]);
    steps     = parse(Int64,   ARGS[3]);
    frec      = parse(Int64,   ARGS[4]);
    param.eta = parse(Float64, ARGS[5]);
    rho       = parse(Float64, ARGS[6]);
    reps      = parse(Float64, ARGS[7]);

elseif length(ARGS) == 0

    # Valores por defecto de parametros

    println("Usando valores por defecto");

    param.N   = 1000
    m_frac    = 0.0006
    steps     = 1000
    frec      = 50
    param.eta = 0.15
    rho       = 10.0
    reps      = 5

elseif length(ARGS) > 0 && length(ARGS) < 7
    println("Faltan parametros")
    exit()
end

tot_links = param.N * (param.N - 1)
links     = m_frac * float(tot_links)

param.M = round(Int, links)
param.bound = sqrt( float(param.N) / rho)

#////////////////////////////////////////////#
#           INITIALIZE SYSTEM                #
#////////////////////////////////////////////#

flock = Flock(param.N, param.M)

T = 1000

out_pos = zeros(Float64, 2 * param.N, T)

# posiciones iniciales
# particulas se inicializan en una caja de tamanio w
flock.pos = Float64[ -param.bound + rand() * 2.0 * param.bound for i in 1:(2 * param.N) ]

# velocidades
# velocidad inicial con direccion aleatoria y norma v0
@time init_vels(param.N, param.v0, flock.vel)

#inicializa red de interaccion no local
nij, poski  = make_IN(param.N, param.M)
flock.nij   = convert(SharedArray, nij)
flock.poski = convert(SharedArray, poski)

ranges     = calc_ranges(param.N)
rij_ranges = calc_ranges(size(flock.rij_ids, 1))

assign_ids(flock.rij_ids, param.N)

flock.vel.s
flock.pos.s
flock.v_r.s
flock.v_n.s
flock.nij.s
flock.poski.s
flock.rij_ids.s
flock.rij.s

param.omega = 1.0
param.eta = 0.0

@time dists(flock, param)
@time dists_id(flock, param, rij_ranges)

p = 2
@profile remotecall(p, loc_vels_chunk_rij, flock.vel, flock.pos, flock.v_r, param.N, param.r0, ranges[p-1])

@sync begin
    for p in procs(flock.pos)
        # @async remotecall(p, loc_vels_chunk, flock.vel, flock.v_r, flock.rij, param.N, range[p-1])
        @async remotecall(p, loc_vels_chunk_rij, flock.vel, flock.pos, flock.v_r, param.N, param.r0, range[p-1])
        @async remotecall(p, non_loc_vels_chunk, flock.vel, flock.v_n, flock.nij, flock.poski, range[p-1])
    end
end

@sync begin
    for p in procs(flock.pos)
        @async remotecall(p, rot_move_chunk, flock.pos, flock.vel, flock.v_r, flock.v_n, flock.noise, dt, param.omega, param.eta, range[p-1])
    end
end


# @time evol_bound(flock, param, out_pos, 100, dt)
@time evol(flock, param, out_pos, T, dt)
@time evol_range(flock, param, out_pos, T, dt, ranges, rij_ranges)

@profile evol_range(flock, param, out_pos, 1, dt, ranges, rij_ranges)

ProfileView.view()

typeof(ranges[1])

print_range(flock.pos)
print_rij_range(param.N)

out_pos

for i in 1:2:param.N
    plt[:plot]( reshape(out_pos[i,:], T), reshape(out_pos[i+1,:], T), "-" )
end

plt[:clf]()

###========================================####

times = [(1 + 10^i):(10^(i+1)) for i in 0:6]
times[1] = 1:10

times

times[1][1] - 1

frecs = [div(length(time), frec) for time in times ]
frecs = [div(length(times[i]), times[i][1]-1) for i in 2:length(times) ]

for k in 1:length(times)
    for t in times[k]

        frec = times[k][1] - 1

        if frec > 0

            if t % frec == 0
                println(t)
            end

        else

            if t % 1 == 0
                println(t)
            end
        end

    end
end
