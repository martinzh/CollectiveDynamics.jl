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
  "-T"
      help     = "Total iterations"
      arg_type = Int64
      default  = 100
  "-s"
      help     = "Sampling Frequency"
      arg_type = Int64
      default  = 10
  "-p"
      help     = "Density -> Particles"
      arg_type = Float64
      default  = 5.0
  "-k"
      help     = "Interaction Network's Connectivity"
      arg_type = Int64
      default  = 2
  "--eta"
      help     = "Noise Intensity"
      arg_type = Float64
      default  = 0.005
  "-w"
      help     = "Relative Weight"
      arg_type = Float64
      default  = 0.005
end

params = parse_args(s)

const dt   = 1.0
const v0   = 1.0

const l    = 0.1 # Regimen de Velocidad

r0 = v0 * dt / l

const L = r0 # Tamaño caja inicial
const N = convert(Int64, L*L* params["p"]) # Numero de particulas (entero)

# const tau      = convert(Int64,2*params["T"])
# const tau      = 750
# const turnRate = 10.0 # velocidad de rotacion
# const ang      = (0.5*pi)/turnRate
# const partPert = RandInt(1,N) # indice de particula perturbada

tau      = convert(Int64,0.5*params["T"])
turnRate = 50.0 # velocidad de rotacion
ang      = (0.5*pi)/turnRate
partPert = RandInt(1,N) # indice de particula perturbada

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

parts    = Array(Bird,N)
dists    = zeros(Float64,N,N) #Matriz de distancias
velPert  = zeros(Float64,2) # vector para velocidad perturbada

println(partPert)
println(tau)

if params["r"] == false
  InitParts(parts,N,10*L,v0,params["k"])
  PrintIntNet(path,parts)
else
  InitParts(parts,N,10*L,v0,params["k"])
  # InitParts(path,N,parts)
  PrintIntNet(path*"_rep",parts)
  println("Passed Replica InitParts")
end

# ==================================== Simulacion ============================================

if params["r"] == false
  for i = 1:params["T"]
    # @time Evoluciona(i,params["step"],parts,params["eta"],params["w"],dists)
    Evoluciona(i,params["s"],parts,params["eta"],params["w"],dists)
  end
else
  for i = 1:params["T"]
    # @time Evoluciona(i,params["step"],parts,params["eta"],params["w"],dists)
    Evoluciona(i,params["s"],parts,params["eta"],params["w"],dists,partPert,velPert,tau,turnRate)
  end
end

# ==================================== Cierra ============================================
close(trays)
close(vels)
