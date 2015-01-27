#################################
# Martin Zumaya Hernandez
# Enero 2015
# Implementacion Flocks en Julia
# Aproximacion de objetos
#################################

# Definicion del tipo Bird

type Bird
  pos::Array{Float64,2}
  vel::Array{Float64,2}
end

########################################################

# Regresa Vector aleatorio entre -L:L

function RandVec(L::Float64)
  x = -L + rand()*2L
  y = -L + rand()*2L
  return [x y]
end

function RandNum(L::Float64)
  x = -L + rand()*2L
  return x
end

########################################################

# Inicializa posiciones y velocidades
# Velocidades normallizadas a v0

function InitParts()
	for i = 1:N
	vel = RandVec(1.0)
	parts[i] = Bird(RandVec(5.0) , scale!(vel,v0/norm(vel)))
	end
end

########################################################

# Actualiza posiciones

function UpdatePos(parts::Array{Bird},dt::Float64)
		for bird = parts
    	# broadcast!(+,bird.pos,bird.pos,scale(bird.vel,dt))
    	bird.pos = bird.pos .+ scale(bird.vel,dt)
		end
end

########################################################

#Construye red aleatoria de conectividad k

function SetLR(k::Int64,X::SparseMatrixCSC{Float64,Int64}) #con sparse
    #Matriz aleatoria
    for i = 1:size(X,1)
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

########################################################

#calcula distancias y adjacencia local

function SetDistAdj(N::Int64,Dist::Array{Float64,2})

  I = Int64[]
  J = Int64[]

	for i = 1:N , j = i+1:N

        d = norm(parts[i].pos-parts[j].pos)

        Dist[i;j] = Dist[j;i] = d

        if d < r0
            push!(I,i)
            push!(J,j)
        end
    end

    return [I J]
end


function SetSR(r0::Float64,Dist::Array{Float64,2},parts::Array{Bird})

    N = size(parts,1)

    I = Int64[]
    J = Int64[]

    for i = 1:N , j = i+1:N

        d = norm(parts[i].pos - parts[j].pos)

        Dist[i;j] = Dist[j;i] = d
        # Dist[i;j] = d

        if d > 0.0 && d < r0
            push!(I,i)
            push!(J,j)
        end
    end
    return sparse(vcat(I,J),vcat(J,I),ones(2*size(I,1)))
end

########################################################

#Calcula las direcciones de la velocidad promedio de cada vecindad
# A -> Matriz de Adjacencia

function GetAngs(parts::Array{Bird}, A::SparseMatrixCSC{Float64,Int64})

	N = size(parts,1)

	# Arreglo de angulos, 1 por particula
    angs = zeros(N)

    # println(angs)

    neigh = rowvals(A) 
    # m, n = size(A)

    # for i = 1:n #itera sobre las particulas con vecindad
    for i = 1:size(A,2) #itera sobre las particulas con vecindad

    		k = 0 # para guardar numero de parts en la vecindad
    		
    		v_prom = parts[i].vel

    		for j = nzrange(A,i)

    			# tal = neigh[j]

    			# println("part : $i ; vecina : $tal")

	        v_prom += parts[neigh[j]].vel 
    			
    			k += 1

    		end

    		# println("k = $k")

        if k > 0
        	scale!(v_prom,1/k)
        	angs[i] = atan2(v_prom[2],v_prom[1]) #agrega el angulo al arreglo
       	end
                                   
    end

    return angs
end

########################################################

#Rota vector 2D

function RotVec(vec::Array{Float64,2},alpha::Float64)

    # M = [cos(alpha) sin(alpha) ; -sin(alpha) cos(alpha)]
    M = [cos(alpha) -sin(alpha) ; sin(alpha) cos(alpha)]

    # return M*vec'
    return vec*M 
end

########################################################

#Actualiza velocidades

function UpdateVel(parts::Array{Bird},SR::SparseMatrixCSC{Float64,Int64},LR::SparseMatrixCSC{Float64,Int64})

    AS = GetAngs(parts,SR) #Angulos inter corto
    AL = GetAngs(parts,LR) #Angulos inter largo

    N = size(parts,1)

    for i = 1:N

    		# dos ruidos, dos intensidades:
        # eta = RandVec(1pi) # Vector de angs aleatorios

        # ang = ruido.*eta # aleatorios * eta

        # ang_tot =  w * (AS[i]+ang[1]) + (1-w) * (AL[i]+ang[2]) 

##############################

        # un ruido, una intensidad

		eta = RandNum(1pi) # angulo aleatorio        

		ang_tot =  w * (AS[i]) + (1-w) * (AL[i]) + eta*ruido

        parts[i].vel = RotVec(parts[i].vel,ang_tot)

    end
end

########################################################

#Actualiza todo

# function Evoluciona(i::Int64, step::Int64, parts::Array{Bird})
function Evoluciona(i::Int64, step::Int64, parts::Array{Bird})

  # @time SR = SetSR(r0,Dist,parts)
  SR = SetSR(r0,Dist,parts)

  UpdatePos(parts,dt)
  UpdateVel(parts,SR,LR)

  # println(i)

  if  i == 1 || i%step == 0
    println("t = $i writing")
    # println(i)
    PrintTrays(i,parts)
    PrintVels(i,parts)
    PrintDist(i,Dist)
  end

end

########################################################
