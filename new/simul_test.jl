# Modelo movimiento colectivo
# Paralelo con SharedArrays

# //////////////// #
# // PARAMETERS // #
# //////////////// #

using CollectiveDynamics

addprocs(4)

N = steps = frec = i = reps = t = 0
rho = m_frac = w = 0.0
dt = v0 = 1.0
l = 0.1

param = Param()
fieldnames(param)

param.r0 = v0 * dt / l
param.w = 0.5; # peso interacciones

if length(ARGS) == 7

    #  Parametros de linea de comando

    N         = parse(Int64,   ARGS[1]);
    m_frac    = parse(Float64, ARGS[2]);
    steps     = parse(Int64,   ARGS[3]);
    frec      = parse(Int64,   ARGS[4]);
    param.eta = parse(Float64, ARGS[5]);
    rho       = parse(Float64, ARGS[6]);
    reps      = parse(Float64, ARGS[7]);

elseif length(ARGS) == 0

    # Valores por defecto de parametros

    println("Usando valores por defecto");

    N         = 1000
    m_frac    = 0.1
    steps     = 1000
    frec      = 50
    param.eta = 0.05
    rho       = 0.1
    reps      = 5

elseif length(ARGS) > 0 && length(ARGS) < 7
    println("Faltan parametros")
    exit()
end

tot_links = N * (N-1)
links     = m_frac * float(tot_links)

param.m = round(Int, links)

w = sqrt( float(N) / rho);
param.bound = w;

#////////////////////////////////////////////#

println(param)

sh_pos  = SharedArray(Float64, 2*N)
sh_vels = SharedArray(Float64, 2*N)

out_pos  = zeros(Float64, 2*N, div(steps, frec) + 1 )
out_vels = zeros(Float64, 2*N, div(steps, frec) + 1 )

# posiciones iniciales
# particulas se inicializan en una caja de tamanio w
sh_pos[:] = Float64[-w + rand() * 2.0 *w for i in 1:2N]

# velocidades
# velocidad inicial con direccion aleatoria y norma v0
for i in 1:2:2N
    vx = -1.0 + 2.0 * rand()
    vy = -1.0 + 2.0 * rand()

    v = norm(Float64[vx,vy])

    sh_vels[i]     = vx / v
    sh_vels[i + 1] = vy / v
end

sh_vels.s
sh_pos.s
