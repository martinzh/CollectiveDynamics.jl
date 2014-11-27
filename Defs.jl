

#Genera N vectores aleatorios de dimension dim
function StartPos()
    for i = 1:N
        #push!(pos,rand(-L:L,dim)) #inicia en una caja de tamanio L
        push!(pos,rand(dim)*L) #inicia en una caja de tamanio L
        push!(vel,rand(dim))
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
function SetMatrix(r0) 

    Adj = zeros(N,N) #Limpia matriz 
    Dist = zeros(N,N) #Limpia matriz 

    for i = 1:N
        for j = N:-1:i
            
            d = norm(pos[i]-pos[j])

            Dist[i;j] = Dist[j;i] = d
            
            #d < r0 && d > 0 ? Adj[i;j] = Adj[j;i] = 1 : Adj[i;j] = 0
            d < r0 ? Adj[i;j] = Adj[j;i] = 1 : Adj[i;j] = 0
        end
    end
end

#Actualiza posiciones
function UpdatePos()

    for i = 1:N
        pos[i] = pos[i] + vel[i] 
    end

end

#Calcula las direcciones de la velocidad promedio de cada vecindad
function GetAngs(A,vel)

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
end

#Construye red aleatoria de conectividad k

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

dim= 2 

dt = 1
f  = 0.1
v0 = 1
ht = 0.3
hg = 0.3
p = 1
l = 0.4
L = 30 # Tamaño caja inicial

N = 4

pos = Array[] #Vector de posiciones
vel = Array[] #Vector de velocidades

Dist = zeros(N,N) #Matriz de distancias

SR_I = zeros(N,N) #Interacciones de corto alcanze
LR_I = zeros(N,N) #Interacciones de lanrgo alcanze -> NO CAMBIA en tiempo

r0 = 3.
L = 2.


StartPos()
StartVels(v0)
SetMatrix(r0)
println("posiciones:")
println(pos)
println("velocidades:")
println(vel)
