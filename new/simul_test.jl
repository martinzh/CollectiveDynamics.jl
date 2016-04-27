# Modelo movimiento colectivo
# Paralelo con SharedArrays

# using CollectiveDynamics
addprocs(4)

@everywhere include("/Users/martinC3/GitRepos/CollectiveDynamics.jl/new/pl.jl")

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

    param.N   = 3500
    m_frac    = 0.1
    steps     = 1000
    frec      = 50
    param.eta = 0.05
    rho       = 0.25
    reps      = 5

elseif length(ARGS) > 0 && length(ARGS) < 7
    println("Faltan parametros")
    exit()
end

tot_links = param.N * (param.N - 1)
links     = m_frac * float(tot_links)

param.M = round(Int, links)
param.bound = sqrt( float(param.N) / rho)

println(param)

#////////////////////////////////////////////#
#           INITIALIZE SYSTEM                #
#////////////////////////////////////////////#

flock = Flock(param.N, param.M)

out_pos = zeros(Float64, 2 * param.N, 1000)

# posiciones iniciales
# particulas se inicializan en una caja de tamanio w
flock.pos = Float64[ -param.bound + rand() * 2.0 * param.bound for i in 1:(2 * param.N) ]

# velocidades
# velocidad inicial con direccion aleatoria y norma v0
init_vels(param.N, param.v0, flock.vel)

flock.vel.s
flock.pos.s

param.omega = 1.0
param.eta = 0.0

# @time evol_bound(flock, param, out_pos, 100, dt)
@time evol(flock, param, out_pos, 1000, dt)

print_range(flock.pos)
print_rij_range(param.N)

out_pos

for i in 1:2:param.N
    plt[:plot]( reshape(out_pos[i,:], 1000), reshape(out_pos[i+1,:], 1000), "-" )
end

plt[:clf]()

###========================================####
