## =========================== ##

# Martin Zumaya Hernandez
# Diciembre 2014
# Simulacion Flocks

## =========================== ##

# Parametros de linea de comandos:

# 1 -> Conectividad
# 2 -> Total de iteraciones
# 3 -> Frecuencia de muestreo
# 4 -> Intensidad de ruido
# 5 -> Peso vecindad geometrica (0 => solo IN ; 1 => solo Geometricas)

## =========================== ##

include("obj_lib.jl")
include("printFunc.jl")

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

if size(ARGS,1) != 0

    const k    = int(ARGS[1])
    const T    = int(ARGS[2])
    const step = int(ARGS[3])
    const eta  = float(ARGS[4])
    const w    = float(ARGS[5])
    const p    = float(ARGS[6])  # Densidad

else
    const k    = 2
    const T    = 500 #iteraciones
    const step = 50 #se recupera informacion cada step
    const eta  = 0.25 #Parametro de ruido
    const w    = 1.0 # Peso relativo de vecindades : 1 => solo IN ; 0 => solo Geometricas
    const p    = 5  # Densidad

end

#Valores default de parametros

const dt   = 1.0
const v0   = 1.0

const l    = 0.1 # Regimen de Velocidad

# const N = convert(Int64, L * L * p) # Numero de particulas (entero)
r0 = v0 * dt / l

const L = r0 # Tamaño caja inicial
const N = int(L * L * p) # Numero de particulas (entero)

# ==================================== Salida de Datos ===========================================

# path = "/Users/martinzh/DATOS_SIMS/DatJul/data_eta$(eta)_k$(k)_w$(w)"
path = "/home/martin/DATOS_SIMS/DatJul/data_eta$(eta)_k$(k)_w$(w)"

MakeDir()
PrintParams()

trays = open("$path/trays.txt","w")
vels = open("$path/vels.txt","w")

println("Particulas = $N")
println("Radio = $r0")
println("Conectividad = $k")

# ==================================== Inicializacion ============================================

parts = Array(Bird,N)
Dist = zeros(N,N) #Matriz de distancias

InitParts(N,10*L,v0,k)

# ==================================== Simulacion ============================================

for i = 1:T
    # @time Evoluciona(i,step,parts,eta,w)
    Evoluciona(i,step,parts,eta,w)
end

# ==================================== Cierra ============================================
close(trays)
close(vels)
