########################

# Martin Zumaya Hernandez
# Diciembre 2014
# Simulacion Flocks

########################


include("flock_lib.jl")

#Escribe trayectorias
function PrintTrays(pos::Array{Array{Float64,1},1})
    for i = 1:N
        @inbounds rr = repr(pos[i])
        write(trays,rr[2:end-1])
        write(trays,"\t")
    end
        write(trays,"\n")
end

#Escribe Matriz Distancias
function PrintDist(i::Int64,Dist::Array{Float64,2})
    d = open("../$path/dists/$i.txt","w")
    for j = 1:size(Dist)[1]
        @inbounds write(d,repr(Dist[j;:])[2:end-1])
        write(d,"\n")
    end
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

    const f = float(ARGS[1])

else
    const f   = 0.099

end

#Valores default de parametros
const dim  = 2 
const dt   = 1.0

const v0   = 1.0
const w    = 0.15 # Peso relativo de vecindades

const hl   = 0.25
const hs   = 0.25

const p    = 1.0
const L    = 25.0 # Tamaño caja inicial
const l    = 0.25

# const T    = 25000 #iteraciones
const T    = 20000 #iteraciones
const step = 400 #se recupera informacion cada step 

const N = convert(Int64, L * L * p) # Numero de particulas (entero)

r0 = v0 * dt / l

ruido = [hl,hs] # [largo,corto]

# println(ruido)

#Crea estructura de folders

path = "../DATA/data_f$(f)"

MakeDir()
PrintParams()

k = convert(Int64,floor(f*N)) # Conectividad en funcion de la fraccion de particulas

trays = open("../$path/trays.txt","w")

println("Particulas = $N")
println("Radio = $r0")
println("Conectividad = $f")

pos = Array{Float64,1}[] #Vector de posiciones
vel = Array{Float64,1}[] #Vector de velocidades

# println(typeof(vel))

Dist = zeros(N,N) #Matriz de distancias

#Usando sparse
LR = spzeros(N,N) #Interacciones de largo alcanze 
                  #No cambia en el tiempo

StartVecs(L,vel,pos)
StartVels(v0,vel)

SetLR(k,LR)

# SR = SetSR(r0,Dist)
# @time GetAngs(vel,SR)
# @time GetAngs(vel,LR)

# @time UpdateVel(vel,SR,LR)
# println(issparse(SR))

for i = 1:T
 
    SR = SetSR(r0,Dist)

    UpdatePos(pos,dt)
    UpdateVel(vel,SR,LR)

    println(i)
    
    if  i == 1 || i%step == 0
        println("writing")
        println(i)
        PrintTrays(pos)
        PrintDist(i,Dist)
    end

end

close(trays)
