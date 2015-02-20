## =========================== ##

# Martin Zumaya Hernandez
# Diciembre 2014
# Simulacion Flocks

## =========================== ##

# Parametros de linea de comandos:

# 1 -> Conectividad
# 2 -> Total de iteraciones
# 3 -> Frecuencia de muestreo

## =========================== ##

include("obj_lib.jl")

# ==================================== Funciones  ==============================================

## =========================== ## ## =========================== ##

#Escribe trayectorias
function PrintTrays(i::Int64,parts::Array{Bird,1})

    # write(trays,"$i\t")

    N = size(parts,1)

    line = ""

    for i in 1:N
        line *= replace(repr(parts[i].pos)[2:end-1],",",'\t')
        if i<N
            line *= "\t"
        end
    end

    line *= "\n"

    write(trays,line)

end

## =========================== ## ## =========================== ##

#Escribe velocidades
function PrintVels(i::Int64,parts::Array{Bird,1})
    
    # write(vels,"$i\t")

    N = size(parts,1)

    line = ""

    for i in 1:N
        line *= replace(repr(parts[i].vel)[2:end-1],",",'\t')
        if i<N
            line *= "\t"
        end
    end

    line *= "\n"

    write(vels,line)

end

## =========================== ## ## =========================== ##

#Escribe Matriz Distancias
function PrintDist(i::Int64,Dist::Array{Float64,2})
    # d = open("../$path/dists/$i.txt","w")
    # d = open("$path/dists/$i.txt","w")
    # writedlm(d,Dist,'\t')
    # close(d)
    writedlm("$path/dists/$i.txt",Dist,'\t')
end

## =========================== ## ## =========================== ##

#Escribe Archivo de Parametros
function PrintParams()
    # d = open("../$path/params.txt","w")
    d = open("$path/params.txt","w")

    write(d,"Particulas = $N\n")
    write(d,"densidad = $p\n")
    write(d,"radio = $r0\n");
    write(d,"f = $f\n");
    write(d,"intensidad de ruido = $eta\n");
    write(d,"peso relativo = $w\n");
    write(d,"regimen de velocidad = $l\n");
    write(d,"iteraciones = $T\n");
    write(d,"frec muestreo = $step\n");
    write(d,"v0 = $v0\n");

    close(d)
end

## =========================== ## ## =========================== ##

function MakeDir()
    try
        # run(`mkdir ../$path`)
        # run(`mkdir ../$path/dists`)

        run(`mkdir $path`)
        run(`mkdir $path/dists`)

    catch y
        println(typeof(y))
    end
end

## =========================== ## ## =========================== ##

# ==================================== Funciones END ============================================

# ==================================== Parametros START ==========================================

# dt  = delta t
# f   = fraccion de conexiones aleatorias relativo al total de parts
# Vo  = magnitud de la velocidad de las particulas
# eta = intensidad del ruido
# pg  = peso de vecindad geometrica
# p   = densidad
# psi = parametro de orden
# l   = regimen de velocidad (tamanio de paso relativo al radio de interaccion)
# r   = radio de interaccion
# N   = numero de particulas

if size(ARGS)[1] != 0

    const f    = float(ARGS[1])
    const T    = int(ARGS[2])
    const step = int(ARGS[3])
    const eta  = float(ARGS[4])

else
    const f    = 0.099
    const T    = 500 #iteraciones
    const step = 50 #se recupera informacion cada step
    const eta  = 0.1 #Parametro de ruido

end

#Valores default de parametros

const dt   = 1.0

const v0   = 1.0
const w    = 0.0 # Peso relativo de vecindades


const p    = 350.0  # Densidad
const L    = 1.0 # Tamaño caja inicial
const l    = 0.25 # Regimen de Velocidad

const N = convert(Int64, L * L * p) # Numero de particulas (entero)

r0 = v0 * dt / l

# ==================================== Parametros END ============================================

# ==================================== Salida de Datos ===========================================

# path = "../DATA/data_f$(f)"

# path = "/home/martin/DATOS_SIMS/DataJul/data_f$(f)"
# path = "/home/martin/DATOS_SIMS/DataJul/data_eta$(eta)"

# path = "/Users/martinzh/DATOS_SIMS/DatJul/data_f$(f)"
path = "/Users/martinzh/DATOS_SIMS/DatJul/data_eta$(eta)"

MakeDir()
PrintParams()

k = convert(Int64,floor(f*N)) # Conectividad en funcion de la fraccion de particulas

trays = open("$path/trays.txt","w")
vels = open("$path/vels.txt","w")

println("Particulas = $N")
println("Radio = $r0")
println("Conectividad = $f")

# ==================================== Inicializacion ============================================

parts = Array(Bird,N)
Dist = zeros(N,N) #Matriz de distancias

#Usando sparse
LR = spzeros(N,N) #Interacciones de largo alcanze
                  #No cambia en el tiempo

InitParts(N,L,v0,k)

# ==================================== Simulacion ============================================

for i = 1:T

    # @time Evoluciona(i,step,parts,eta,w)
    Evoluciona(i,step,parts,eta,w)

end

# ==================================== Cierra ============================================

close(trays)
close(vels)
