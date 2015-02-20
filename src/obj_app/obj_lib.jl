## =========================== ##

# Martin Zumaya Hernandez
# Enero 2015
# Implementacion Flocks en Julia
# Aproximacion de objetos

## =========================== ##

# Definicion del tipo Bird

# type Bird
#   pos::Array{Float64,1}
#   vel::Array{Float64,1}
# end

type Bird
  pos::Array{Float64,1}
  vel::Array{Float64,1}
  inputs::Array{Int64,1}
end

## =========================== ## ## =========================== ##

# Regresa Vector aleatorio entre -L:L

function RandVec(L::Float64)
  return [-L + rand()*2.0*L, -L + rand()*2.0*L]
end

function RandNum(L::Float64)
  return -L + rand()*2.0*L
end

## =========================== ## ## =========================== ##

# Inicializa posiciones y velocidades
# Velocidades normallizadas a v0

# function InitParts(N::Int64,L::Float64,v0::Float64)
#     for i = 1:N
#     vel = RandVec(1.0)
#     parts[i] = Bird(RandVec(L) , scale!(vel,v0/norm(vel)))
#     end
# end

function InitParts(N::Int64,L::Float64,v0::Float64,k::Int64)
	for i = 1:N
	vel = RandVec(1.0)
	parts[i] = Bird(RandVec(L) , scale!(vel,v0/norm(vel)), Inputs(k,N))
	end
end

## =========================== ## ## =========================== ##

# Actualiza posiciones

function UpdatePos!(parts::Array{Bird,1},dt::Float64)
		for bird = parts
    	broadcast!(+,bird.pos,bird.pos,scale(bird.vel,dt))
    	# bird.pos = bird.pos .+ scale(bird.vel,dt)
		end
end

## =========================== ## ## =========================== ##

#Construye red aleatoria de conectividad k

function Inputs(k::Int64,N::Int64)
    
    ins = Array(Int64,k)
    
    for j = 1:k
        switch = true
            while switch
                s = rand(1:N)
                if s != j
                    ins[j] = s
                    switch = false
                end
            end
    end 
    return ins
end

## =========================== ## ## =========================== ##

function SetLR(k::Int64,N::Int64,X::SparseMatrixCSC{Float64,Int64}) 

    for i = 1:size(X,1) , j = 1:k
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

## =========================== ## ## =========================== ##

#calcula distancias y adjacencia local

function SetSR(r0::Float64,Dist::Array{Float64,2},parts::Array{Bird,1})

    N = size(parts,1)

    I = Int64[]
    J = Int64[]

    for i = 1:N , j = i+1:N

        d = norm(parts[i].pos - parts[j].pos)

        Dist[i,j] = Dist[j,i] = d
        # Dist[i;j] = d

        # if d > 0.0 && d < r0
        if d < r0
            push!(I,i)
            push!(J,j)
        end
    end
    return sparse(vcat(I,J),vcat(J,I),ones(2*size(I,1)))
end

## =========================== ## ## =========================== ##

#Calcula las direcciones de la velocidad promedio de cada vecindad
# A -> Matriz de Adjacencia

function GetAngs(parts::Array{Bird,1}, A::SparseMatrixCSC{Float64,Int64})

    N = size(parts,1)

    # Arreglo de angulos, 1 por particula
    angs = zeros(N)

    # println(angs)

    neigh = rowvals(A)

    # for i = 1:n #itera sobre las particulas con vecindad
    for i = 1:size(A,2) #itera sobre las particulas con vecindad

            k = 0.0 # para guardar numero de parts en la vecindad

            v_prom = parts[i].vel

            for j = nzrange(A,i)

                # tal = neigh[j]

                # println("part : $i ; vecina : $tal")

               v_prom += parts[neigh[j]].vel

                k += 1.0

            end

            # println("k = $k")

        if k > 0
            scale!(v_prom,1.0/k)
            angs[i] = atan2(v_prom[2],v_prom[1]) #agrega el angulo al arreglo
        end

    end

    return angs
end

## =========================== ## ## =========================== ##

# Angulos de entradas en la red de itneraccion

function GetAngsIN(parts::Array{Bird,1})

	k = size(parts[1].inputs,1) # conectividad

    N = size(parts,1)

	# Arreglo de angulos, 1 por particula
    angs = zeros(N)

    # println(angs)

    if k>0

        for i = 1:N
            
            v_prom = parts[i].vel

            for j = parts[i].inputs
                v_prom += parts[j].vel
            end

            scale!(v_prom,1.0/k)
            angs[i] = atan2(v_prom[2],v_prom[1]) #agrega el angulo al arreglo
        end

    end

    return angs
end

## =========================== ## ## =========================== ##

#Rota vector 2D

function RotVec!(vec::Array{Float64,1},alpha::Float64)

    X = vec[1]*cos(alpha) - vec[2]*sin(alpha)
    Y = vec[1]*sin(alpha) + vec[2]*cos(alpha)

    vec[1] = X
    vec[2] = Y
end

## =========================== ## ## =========================== ##

#Actualiza velocidades

function UpdateVel!(parts::Array{Bird,1},SR::SparseMatrixCSC{Float64,Int64},LR::SparseMatrixCSC{Float64,Int64},ruido::Float64,w::Float64)

    AS = GetAngs(parts,SR) #Angulos inter corto
    AL = GetAngsIN(parts) #Angulos inter largo

    for i = 1:size(parts,1)

      ang_rand = RandNum(1.0*pi) # angulo aleatorio

      ang_tot =  w * (AS[i]) + (1-w) * (AL[i]) + ang_rand*ruido

      RotVec!(parts[i].vel,ang_tot)

    end
end

## =========================== ## ## =========================== ##

#Actualiza todo

function Evoluciona(i::Int64, step::Int64, parts::Array{Bird,1},ruido::Float64,w::Float64)

  # @time SR = SetSR(r0,Dist,parts)
  SR = SetSR(r0,Dist,parts)

  UpdatePos!(parts,dt)
  UpdateVel!(parts,SR,LR,ruido,w)

  # println(i)

  if  i == 1 || i%step == 0
    println("t = $i writing")
    # println(i)
    PrintTrays(i,parts)
    PrintVels(i,parts)
    PrintDist(i,Dist)
  end

end

## =========================== ## ## =========================== ##
