# Modelo movimiento colectivo
# Paralelo con SharedArrays

# using CollectiveDynamics

@everywhere reload("CollectiveDynamics.jl")

addprocs(4)

using PyPlot

####========================================####
####      SHARED ARRAYS UILITY FUNCTIONS    ####
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
    splits[idx]+1:splits[idx+1]
end

####===============================####

@everywhere function rot_move_chunk(pos::SharedArray, vel::SharedArray, id_range::UnitRange, dt::Float64)

    for i in id_range
        pos[i] += vel[i] * dt
    end

end

@everywhere rot_move(pos::SharedArray, vel::SharedArray, dt::Float64) = rot_move_chunk(pos, vel, myrange(pos), dt)
####===============================####

function evol(pos::SharedArray, vel::SharedArray, vals::Array{Float64,2}, T::Int64, dt::Float64)
    for t in 1:T
        @sync begin
            for p in procs(pos)
                @async remotecall_wait(p, rot_move, pos, vel, dt)
            end
        end
        vals[:,t] = pos.s
    end
end
####===============================####

# function advection_shared!(q, u)
#     @sync begin
#         for p in procs(q)
#             @async remotecall_wait(p, advection_shared_chunk!, q, u)
#         end
#     end
#     q
# end

function print_range(q)
    @sync begin
        for p in procs(q)
            # @async remotecall_wait(p, advection_shared_chunk!, q, u)
            @async remotecall_wait(p, myrange, q)
        end
    end

end

####========================================####
####             PARAMETERS                 ####
####========================================####

N = steps = frec = i = reps = t = 0
rho = m_frac = w = 0.0
dt = v0 = 1.0
# const dt = 1.0
l = 0.1

param = CollectiveDynamics.Param()
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

w = sqrt( float(N) / rho)
param.bound = w

#////////////////////////////////////////////#

println(param)

sh_pos     = SharedArray(Float64, 2N) # posiciones
sh_vels    = SharedArray(Float64, 2N) # velocidades
sh_loc     = SharedArray(Float64, N) # angulos locales
sh_non_loc = SharedArray(Float64, N) # angulos no locales

sh_rij = convert(SharedArray, zeros(N,N))

# out_pos  = zeros(Float64, 2*N, div(steps, frec) + 1 )
# out_vels = zeros(Float64, 2*N, div(steps, frec) + 1 )

out_pos = zeros(Float64, 2N, 10)

# posiciones iniciales
# particulas se inicializan en una caja de tamanio w
sh_pos[:] = Float64[-w + rand() * 2.0 *w for i in 1:2N]

# velocidades
# velocidad inicial con direccion aleatoria y norma v0
CollectiveDynamics.init_vels(N, v0, sh_vels)

sh_vels.s
sh_pos.s

print_range(sh_pos)

@time evol(sh_pos, sh_vels, out_pos, 10, dt)

out_pos

# plt[:plot](Float64[out_pos[i,1] for i in 1:2:2N], Float64[out_pos[i,1] for i in 2:2:2N], ".-" )

for i in 1:2:N
    plt[:plot]( reshape(out_pos[i,:], 10), reshape(out_pos[i+1,:],10), "-" )
end

plt[:clf]()

for i in 2:2:10
    println(i)
end

Float64[out_pos[i,1] for i in 2:2:2N]


rij = Dict()

for i in 1:5, j in (i+1):5
    # println(i,",",j)
    rij[(i,j)] = rand(0:1)
end

rij

m = 5

i = 3

u = 0

for j in 1:(i-1)
    # println(i,",",j)
    println((j,i),",",rij[(j,i)])
    u += rij[(j,i)]
end

for j in (i+1):m
    # println(i,",",j)
    println((i,j),",",rij[(i,j)])
    u += rij[(i,j)]
end

u
