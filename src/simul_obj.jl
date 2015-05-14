## =========================== ##

# Martin Zumaya Hernandez
# Diciembre 2014
# Simulacion Flocks

## =========================== ##

include("obj_lib.jl")
include("printFunc.jl")

using ArgParse

# ==================================== Parametros START ==========================================

# dt  = delta t
# k   = conectividad de la red de interaccion
# Vo  = magnitud de la velocidad de las particulas
# eta = intensidad del ruido
# pg  = peso de vecindad geometrica
# p   = densidad
# psi = parametro de orden
# l   = regimen de velocidad (tamanio de paso relativo al radio de interaccion)
# r   = radio de interaccion
# N   = numero de particulas

s = ArgParseSettings()

@add_arg_table s begin
  "-m"
      help    = "selects OS"
      action  = :store_true
      default = false
  "-r"
      help    = "replica system"
      action  = :store_true
      default = false
  "k"
      help     = "Interaction Network's Connectivity"
      arg_type = Int64
      default  = 2
  "T"
      help     = "Total iterations"
      arg_type = Int64
      default  = 100
  "step"
      help     = "Sampling Frequency"
      arg_type = Int64
      default  = 10
  "eta"
      help     = "Noise Intensity"
      arg_type = Float64
      default  = 0.005
  "w"
      help     = "Relative Weight"
      arg_type = Float64
      default  = 0.005
  "p"
      help     = "Density -> Particles"
      arg_type = Float64
      default  = 5.0
end

params = parse_args(s)

const dt   = 1.0
const v0   = 1.0

const l    = 0.1 # Regimen de Velocidad

r0 = v0 * dt / l

const L = r0 # Tamaño caja inicial
const N = int(L * L * params["p"]) # Numero de particulas (entero)

# ==================================== Salida de Datos ===========================================

if params["m"] == true
  path = "/Users/martinzh/DATOS_SIMS/DatJul/data_eta$(params["eta"])_k$(params["k"])_w$(params["w"])"
else
  path = "/home/martin/DATOS_SIMS/DatJul/data_eta$(params["eta"])_k$(params["k"])_w$(params["w"])"
end

println("Particulas   = $N")
println("Radio        = $r0")
println("Conectividad = $(params["k"])")


if params["r"] == false

  MakeDir(path)
  PrintParams(path)

  trays = open("$path/trays.txt","w")
  vels  = open("$path/vels.txt","w")

else

  MakeDir(path*"_rep")
  PrintParams(path*"_rep")

  trays = open("$(path)_rep/trays.txt","w")
  vels  = open("$(path)_rep/vels.txt","w")

end

# ==================================== Inicializacion ============================================

parts = Array(Bird,N)
dists = zeros(Float64,N,N) #Matriz de distancias

if params["r"] == false
  InitParts(parts,N,10*L,v0,params["k"])
else
  InitParts(path,parts)
end

PrintIntNet(parts)

# ==================================== Simulacion ============================================

for i = 1:params["T"]
    @time Evoluciona(i,params["step"],parts,params["eta"],params["w"],dists)
    Evoluciona(i,params["step"],parts,params["eta"],params["w"],dists)
end

# ==================================== Cierra ============================================
close(trays)
close(vels)
