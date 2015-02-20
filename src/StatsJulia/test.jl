include("StatsLib.jl")
using Gadfly

eta = 0.0
path = "/Users/martinzh/DATOS_SIMS/DatJul/data_eta$eta"

params = Tmp.GetParams(path)
vels = Tmp.GetVecs(path,"vels")
trays = Tmp.GetVecs(path,"trays")

println(typeof(vels))
println(params)

N = int(params[1])
ro = params[2]
radio = params[3]
f = params[4]
eta = params[5]
relWeight = params[6]
regVel = params[7]
tf = int(params[8])
step = int(params[9])
v0 = params[10]

itTot = int(tf/step)

#vels
#vels[1,1:end]

tiempo = vcat([1],[ i*step for i = 1:itTot])

#typeof(tiempo[1])

velsProm = Tmp.VecProm(vels,itTot,step)
posProm = Tmp.VecProm(trays,itTot,step)

size(velsProm)


ordParam = Tmp.OrderParam(velsProm,v0)

plot(x = tiempo, y = ordParam,
     Geom.point, Geom.line,
     Guide.xlabel("Tiempo"), Guide.ylabel("Ψ"),Guide.title("Parámetro de orden \n η = $eta")
     )
