

#Genera N vectores aleatorios de dimension dim
function StartVecs(L::Float64,vel::Array{Array{Float64,2},1},pos::Array{Array{Float64,2},1})
    for i = 1:N
        push!(pos,[rand(-L:L) rand(-L:L)]) #inicia en una caja de tamanio L
        push!(vel,[rand() rand()])
    end
end

#Normaliza los vectores, magnitud v0 direccion aleatoria
function StartVels(v0::Float64,vel::Array{Array{Float64,2},1})
    for i in 1:N
        @inbounds vel[i] = vel[i] * (v0/norm(vel[i]))
        #println(norm(vel[i]))
    end
end

#calcula distancias y adjacencia local
# function SetMatrix(r0::Float64,SR::Array{Float64,2},Dist::Array{Float64,2})
function SetMatrix(r0::Float64,SR::SparseMatrixCSC{Float64,Int64},Dist::Array{Float64,2})
    for i = 1:N
        for j = N:-1:i
            
            @inbounds d = norm(pos[i]-pos[j])

            @inbounds Dist[i;j] = Dist[j;i] = d
            
            #d < r0 && d > 0 ? SR[i;j] = SR[j;i] = 1 : SR[i;j] = 0
            @inbounds d < r0 ? SR[i;j] = SR[j;i] = 1 : SR[i;j] = 0
        end
    end
end

#Actualiza posiciones
function UpdatePos(pos::Array{Array{Float64,2},1})
    for i = 1:N
        @inbounds pos[i] = pos[i] + vel[i]*dt
    end
end

#Calcula las direcciones de la velocidad promedio de cada vecindad
# function GetAngs(A,vel)
# function GetAngs(A::Array{Float64,2})
function GetAngs(A::SparseMatrixCSC{Float64,Int64}) #con sparse

    angs = Float64[] #Arreglo para guardar angulos

    # println(angs)
    # Calcula los angulos de las velocidades promedio

    for i = 1:size(A)[1]

        # U = vel'.*A[i;:] #Determina la vecindad
        U = sum(broadcast(*,vel',A[i;:])) #com sparse

        k = sum(A[i;:]) #Numero de parts en la vecindad 

        # v_prom = 1/k * sum(U) #Calcula vector promedio
        v_prom = 1/k * U #Calcula vector promedio

        # # println("part $i -> ang prom $(atan2(v_prom[2],v_prom[1]))")

        push!(angs,atan2(v_prom[2],v_prom[1])) #agrega el angulo al arreglo
    end
    # println("Angulos:\n $angs")
    return angs
end

#Construye red aleatoria de conectividad k
# function SetLR(k::Int64,X::Array{Float64,2})
function SetLR(k::Int64,X::SparseMatrixCSC{Float64,Int64}) #con sparse
    #Matriz aleatoria
    for i = 1:size(X)[1]
        for j = 1:k        
            switch = true
            while switch 
                s = rand(1:N)
                if s != i
                    X[i,s] = 1
                    #X[i,s] = X[s,i] = 1 
                    switch = false
                end
            end
        end
    end
end

#Actualiza velocidades
# function UpdateVel(vel::Array{Array{Float64,2},1},SR::Array{Float64,2},LR::Array{Float64,2})
function UpdateVel(vel::Array{Array{Float64,2},1},SR::SparseMatrixCSC{Float64,Int64},LR::SparseMatrixCSC{Float64,Int64})

    AS = GetAngs(SR) #Angulos inter corto
    AL = GetAngs(LR) #Angulos inter largo

    # println("Corto Alcance:")
    # println(AS)
    # println("Largo Alcance:")
    # println(AL)

    for i = 1:N

        eta = rand(-pi:pi,2) # Vector de angs aleatorios

        # println(ruido)
        # println(eta)

        ang = ruido'.*eta # aleaorios * eta
        # println(ang)

        ang_tot = abs(1-w)*(AL[i]+ang[1]) + w*(AS[i]+ang[2])

        vel_n = RotVec(vel[i],ang_tot)
        
        # println("original")
        # println(vel[i])

        # println("nuevo")
        # println(vel_n)

        vel[i] = vel_n

        # println("rotado")
        # println(vel[i])

    end
end

#Rota vector 2D
function RotVec(vec::Array{Float64,2},alpha::Float64)

    # println(vec)

    M = [cos(alpha) sin(alpha) ; -sin(alpha) cos(alpha)]

    res = M*vec'

    # println(typeof(res))

    # println(res)

    return res'
end

#Escribe trayectorias
function PrintTrays(pos::Array{Array{Float64,2},1})
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
    write(d,"ruido geometrico = $hg\n");
    write(d,"ruido topologico = $ht\n");
    write(d,"peso geometrico = $w\n");
    write(d,"regimen de velocidad = $l\n");
    write(d,"iteraciones = $T\n");

    close(d)
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

const ht   = 0.25
const hg   = 0.25
const p    = 1.25
const L    = 20.0 # Tamaño caja inicial
const l    = 0.35
# const T    = 25000 #iteraciones
const T    = 1000 #iteraciones
const step = 250 #se recupera informacion cada step 



const N = convert(Int64, L*L * p) # Numero de particulas (entero)

r0 = v0 * dt / l

ruido = [ht hg] # [largo corto]

# println(ruido)

#Crea estructura de folders

# path = "data_f$(f)_w$w"
# path = "JuliaFlocks/Julia/DATA/data_f$(f)"
path = "../DATA/data_f$(f)"

# run(`if [ ! -d ""../$path""  ]; then
#   mkdir ../$path
#   mkdir ../$path/dists
#     fi`)

run(`mkdir ../$path`)
run(`mkdir ../$path/dists`)
# run(`mkdir ../$path/adjs`)

PrintParams()

k = convert(Int64,floor(f*N)) # Conectividad en funcion de la fraccion de particulas

trays = open("../$path/trays.txt","w")

println("Particulas = $N")
println("Radio = $r0")
println("Conectividad = $f")


pos = Array{Float64,2}[] #Vector de posiciones
vel = Array{Float64,2}[] #Vector de velocidades

# println(typeof(vel))

Dist = zeros(N,N) #Matriz de distancias

# SR = zeros(N,N) #Interacciones de corto alcanze
# LR = zeros(N,N) #Interacciones de lanrgo alcanze -> NO CAMBIA en tiempo

#Usando sparse
SR = spzeros(N,N) #Interacciones de corto alcanze
LR = spzeros(N,N) #Interacciones de lanrgo alcanze -> NO CAMBIA en tiempo

StartVecs(L,vel,pos)
StartVels(v0,vel)

# println(typeof(k))
# println(typeof(LR))
# println(typeof(SR))

SetLR(k,LR)

# println("$vel\n")

# for i = 1:N

#     println("original")
#     println(vel[i])

#     otro = Array(2*vel[i])
    
#     println("otro")
#     println(otro)

#     vel[i] = otro

#     # vel[i] = 2*vel[i]

#     println("asignado")
#     println(vel[i])
# end

for i = 1:T

    # println("posiciones:\n $pos")
    # println("velocidades:\n $vel")

    SetMatrix(r0,SR,Dist)
    UpdatePos(pos)
    UpdateVel(vel,SR,LR)
    
    if  i == 1 || i%step == 0
        println(i)
        PrintTrays(pos)
        PrintDist(i,Dist)
    end

    # println("Long Range:\n $LR")
    # println("Short Range:\n $SR")
    # println("Distancias:\n $Dist")

end

close(trays)
