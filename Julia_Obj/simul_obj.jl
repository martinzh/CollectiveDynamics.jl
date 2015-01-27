########################

# Martin Zumaya Hernandez
# Diciembre 2014
# Simulacion Flocks

########################

########################

# Parametros de linea de comandos:

# 1 -> Conectividad
# 2 -> Total de iteraciones
# 3 -> Frecuencia de muestreo

########################

include("obj_lib.jl")

#Escribe trayectorias
function PrintTrays(i::Int64,parts::Array{Bird})

    write(trays,"$i\t")
    
    for bird = parts
        write(trays,repr(bird.pos)[2:end-1])
        # writedlm(trays,bird.pos,'\t')
        # write(trays,bird.pos)
        write(trays,"\t")
    end

    write(trays,"\n")
end

#Escribe velocidades
function PrintVels(i::Int64,parts::Array{Bird})
    
    write(vels,"$i\t")

    for bird = parts
        write(vels,repr(bird.vel)[2:end-1])
        # write(vels,bird.vel,'\t')
        # write(vels,bird.vel)
        write(vels,"\t")
    end

    write(vels,"\n")
end

#Escribe Matriz Distancias
function PrintDist(i::Int64,Dist::Array{Float64,2})
    d = open("../$path/dists/$i.txt","w")
    # for j = 1:size(Dist)[1]
    #     write(d,repr(Dist[j;:])[2:end-1])
    #     write(d,"\n")
    # end
    writedlm(d,Dist,'\t')
    close(d)
end

#Escribe Archivo de Parametros
function PrintParams()
    d = open("../$path/params.txt","w")

    write(d,"Particulas = $N\n")
    write(d,"densidad = $p\n")
    write(d,"radio = $r0\n");
    write(d,"f = $f\n");
    write(d,"ruido corto alcance = $hs\n");
    write(d,"ruido largo alcance = $hl\n");
    write(d,"peso relativo = $w\n");
    write(d,"regimen de velocidad = $l\n");
    write(d,"iteraciones = $T\n");
    write(d,"frec muestreo = $step\n");

    close(d)
end

function MakeDir()
    try
        run(`mkdir ../$path`)
        run(`mkdir ../$path/dists`)
        # run(`mkdir ../$path/adjs`)
    catch y
        println(typeof(y))
    end
end


# ==================================== Parametros ==============================================

    # dt -> delta t
    # f -> fraccion de conexiones aleatorias relativo al total de parts
    # Vo -> magnitud de la velocidad de las particulas
    # ht ->ruido topologico
    # hg -> ruido geometrico
    # pg -> peso de vecindad geometrica
    # p  -> densidad
    # psi -> parametro de orden
    # l  -> regimen de velocidad (tamanio de paso relativo al radio de interaccion)
    # r -> radio de interaccion
    # N -> numero de particulas

if size(ARGS)[1] != 0

    const f    = float(ARGS[1])
    const T    = int(ARGS[2])
    const step = int(ARGS[3])

else
    const f    = 0.099
    const T    = 25000 #iteraciones
    const step = 500 #se recupera informacion cada step

end

#Valores default de parametros
const dim  = 2
const dt   = 1.0

const v0   = 1.0
const w    = 0.15 # Peso relativo de vecindades

const hl   = 0.1
const hs   = 0.1

const p    = 1.2  # Densidad
const L    = 30.0 # Tamaño caja inicial
const l    = 0.25 # Regimen de Velocidad


const N = convert(Int64, L * L * p) # Numero de particulas (entero)

r0 = v0 * dt / l

# ruido = [hl,hs] # [largo,corto]

ruido = hl # solo un parametro de intensidad de ruido

# println(ruido)

#Crea estructura de folders

path = "../DATA/data_f$(f)"

MakeDir()
PrintParams()

k = convert(Int64,floor(f*N)) # Conectividad en funcion de la fraccion de particulas

trays = open("../$path/trays.txt","w")
vels = open("../$path/vels.txt","w")

println("Particulas = $N")
println("Radio = $r0")
println("Conectividad = $f")

parts = Array(Bird,N)

Dist = zeros(N,N) #Matriz de distancias

#Usando sparse
LR = spzeros(N,N) #Interacciones de largo alcanze
                  #No cambia en el tiempo

InitParts()

SetLR(k,LR)

for i = 1:T

    @time Evoluciona(i,step,parts)
    # Evoluciona(i,step,parts)

end

close(trays)
