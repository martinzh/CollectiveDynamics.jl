

#Genera N vectores aleatorios de dimension dim
function StartPos()
    for i = 1:N
        LL = convert(Float64,L)
        # push!(pos,rand(-LL:LL,dim)) #inicia en una caja de tamanio L
        # push!(vel,rand(dim))
        push!(pos,[rand(-LL:LL) rand(-LL:LL)]) #inicia en una caja de tamanio L
        push!(vel,[rand() rand()])
    end
end

#Normaliza los vectores, magnitud v0 direccion aleatoria
function StartVels(v0)
    for i in 1:N
        vel[i] = vel[i] * (v0/norm(vel[i]))
        #println(norm(vel[i]))
    end
end

#calcula distancias y adjacencia local
function SetMatrix(r0,SR,Dist) 

    # SR = zeros(N,N) #Limpia matriz 
    # Dist = zeros(N,N) #Limpia matriz 

    for i = 1:N
        for j = N:-1:i
            
            d = norm(pos[i]-pos[j])

            Dist[i;j] = Dist[j;i] = d
            
            #d < r0 && d > 0 ? SR[i;j] = SR[j;i] = 1 : SR[i;j] = 0
            d < r0 ? SR[i;j] = SR[j;i] = 1 : SR[i;j] = 0
        end
    end
end

#Actualiza posiciones
function UpdatePos()

    for i = 1:N
        pos[i] = pos[i] + vel[i]*dt
    end

end

#Calcula las direcciones de la velocidad promedio de cada vecindad
# function GetAngs(A,vel)
function GetAngs(A)

    angs = Float64[] #Arreglo para guardar angulos

    # println(angs)
    # Calcula los angulos de las velocidades promedio

    for i = 1:size(A)[1]

        U = vel'.*A[i;:] #Determina la vecindad
        k = sum(A[i;:]) #Numero de parts en la vecindad 

        v_prom = 1/k * sum(U) #Calcula vector promedio
    #    println(degs(atan2(v[2],v[1])))

        # a = atan2(v_prom[2],v_prom[1]) #Direccion del vector promedio
        # push!(angs,a) #agrega el angulo al arreglo
        # # println("part $i -> ang prom $a")

        push!(angs,atan2(v_prom[2],v_prom[1])) #agrega el angulo al arreglo
    end
    # println("Angulos:\n $angs")
    return angs
end

#Construye red aleatoria de conectividad k
function SetLR(k,X)
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
function UpdateVel()

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
function RotVec(vec,alpha)

    # println(vec)

    M = [cos(alpha) sin(alpha) ; -sin(alpha) cos(alpha)]

    res = M*vec'

    # println(typeof(res))

    # println(res)

    return res'
end

#Escribe trayectorias
function PrintTrays()
    for i = 1:N
        rr = repr(pos[i])
        write(trays,rr[2:end-1])
        write(trays,"\t")
    end
        write(trays,"\n")
end

#Escribe Matriz Distancias
function PrintDist(i)
    d = open("../$path/dists/$i.txt","w")
    for j = 1:size(Dist)[1]
        write(d,repr(Dist[j;:])[2:end-1])
        write(d,"\n")
    end
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

#Valores default de parametros

dim = 2 

dt  = 1
f   = 0.05
v0  = 1
ht  = 0.25
hg  = 0.25
p   = 0.1
l   = 0.25
L   = 30 # Tamaño caja inicial
w   = 0.3 # Peso relativo de vecindades

T   = 10 #iteraciones

N = convert(Int64, L^2 * p) # Numero de particulas (entero)

r0 = v0 * dt / l

ruido = [ht hg] # [largo corto]

# println(ruido)

#Crea estructura de folders

path = "data_f$(f)_ro$p"

run(`mkdir ../$path`)
run(`mkdir ../$path/dists`)

f = convert(Int64,floor(f*N)) # Conectividad en funcion de la fraccion de particulas

trays = open("../$path/trays.txt","w")

println("Particulas = $N")
println("Radio = $r0")
println("Conectividad = $f")

pos = Array[] #Vector de posiciones
vel = Array[] #Vector de velocidades

Dist = zeros(N,N) #Matriz de distancias

SR = zeros(N,N) #Interacciones de corto alcanze
LR = zeros(N,N) #Interacciones de lanrgo alcanze -> NO CAMBIA en tiempo

StartPos()
StartVels(v0)
SetLR(f,LR)

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

#     # println("posiciones:\n $pos")
#     # println("velocidades:\n $vel")
    println(i)
    SetMatrix(r0,SR,Dist)
    UpdatePos()
    UpdateVel()
    
    if  i%2 == 0

        PrintTrays()
        PrintDist(i)
    end
#     # println("Long Range:\n $LR")
#     # println("Short Range:\n $SR")
#     # println("Distancias:\n $Dist")

end

close(trays)
