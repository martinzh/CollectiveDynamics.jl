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

if size(ARGS)[1] != 0

    const k    = int(ARGS[1])
    const T    = int(ARGS[2])
    const step = int(ARGS[3])
    const eta  = float(ARGS[4])
    const w    = float(ARGS[5])

else
    const k    = 2
    const T    = 500 #iteraciones
    const step = 50 #se recupera informacion cada step
    const eta  = 0.1 #Parametro de ruido
    const w    = 1.0 # Peso relativo de vecindades : 0 => solo IN ; 1 => solo Geometricas

end

#Valores default de parametros

const dt   = 1.0

const v0   = 1.0

const p    = 5.0  # Densidad
const l    = 0.1 # Regimen de Velocidad

# const N = convert(Int64, L * L * p) # Numero de particulas (entero)
r0 = v0 * dt / l

const L = r0 # Tamaño caja inicial
const N = int(L * L * p) # Numero de particulas (entero)


# ==================================== Parametros END ============================================

# ==================================== Salida de Datos ===========================================

<<<<<<< HEAD
# path = "/home/martin/DATOS_SIMS/DataJul/data_eta$(eta)_k$(k)_w$(w)"
path = "/home/martin/DATOS_SIMS/DataJul/data_eta$(eta)_k$(k)"

# path = "/Users/martinzh/DATOS_SIMS/DatJul/data_eta$(eta)_k$(k)_w$(w)"
# path = "/Users/martinzh/DATOS_SIMS/DatJul/data_eta$(eta)_k$(k)"
=======
# path = "/home/martin/DATOS_SIMS/DataJul/data_eta$(eta)_k$(k)"

path = "/Users/martinzh/DATOS_SIMS/DatJul/data_eta$(eta)_k$(k)"
>>>>>>> d72f798f8b217669cf805ccbbcc5afb160c22f5e

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

InitParts(N,4*L,v0,k)

#Usando sparse ==> Se hace con los arreglos de cada particula
# LR = spzeros(N,N) #Interacciones de largo alcanze
                  #No cambia en el tiempo


# ==================================== Simulacion ============================================

for i = 1:T

    # @time Evoluciona(i,step,parts,eta,w)
    Evoluciona(i,step,parts,eta,w)

end

# ==================================== Cierra ============================================

close(trays)
close(vels)
