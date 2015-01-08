########################

# Martin Zumaya Hernandez
# Diciembre 2014
# Libreria para flocks

########################

#Genera N vectores aleatorios de dimension dim

function StartVecs(L::Float64, vel::Array{Array{Float64,1},1}, pos::Array{Array{Float64,1},1})
    for i = 1:N
        push!(pos,[rand(-L:L),rand(-L:L)]) #inicia en una caja de tamanio L
        push!(vel,[rand(),rand()])
    end
end

########################################################

#Normaliza los vectores, magnitud v0 direccion aleatoria

function StartVels(v0::Float64, vel::Array{Array{Float64,1},1})
    for i in 1:N
        scale!(vel[i],v0/norm(vel[i])) #escala y modifica vel[i]
    end
end

########################################################

#Actualiza posiciones

function UpdatePos(pos::Array{Array{Float64,1},1},dt::Float64)
    broadcast!(+,pos,pos,scale(vel,dt))
end

########################################################

#Construye red aleatoria de conectividad k

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

########################################################

#calcula distancias y adjacencia local
function SetSR(r0::Float64,Dist::Array{Float64,2})

    I = Int64[]
    J = Int64[]

    for i = 1:N , j = N:-1:i

        d = norm(pos[i]-pos[j])

        Dist[i;j] = Dist[j;i] = d
        # Dist[j;i] = d

        if d < r0
            push!(I,i)
            push!(J,j)
        end
    end

    # Dist += Dist'

    return sparse(vcat(I,J),vcat(J,I),ones(2*size(I)[1]))
end

########################################################

#Calcula las direcciones de la velocidad promedio de cada vecindad
function GetAngs(vel::Array{Array{Float64,1},1}, A::SparseMatrixCSC{Float64,Int64}) #con sparse

    angs = Array(Float64,size(A)[1])

    for i = 1:size(A)[1]

        U = sum(broadcast(*,vel,A[:;i])) #com sparse

        k = nnz(A[:;i]) #Numero de parts en la vecindad

        scale!(U,1/k)

        angs[i] = atan2(U[2],U[1]) #agrega el angulo al arreglo
    end

    return angs
end

########################################################

#Rota vector 2D

function RotVec(vec::Array{Float64,1},alpha::Float64)

    M = [cos(alpha) sin(alpha) ; -sin(alpha) cos(alpha)]

    return M*vec
end

########################################################

#Actualiza velocidades

function UpdateVel(vel::Array{Array{Float64,1},1},SR::SparseMatrixCSC{Float64,Int64},LR::SparseMatrixCSC{Float64,Int64})

    AS = GetAngs(vel,SR) #Angulos inter corto
    AL = GetAngs(vel,LR) #Angulos inter largo

    for i = 1:N

        eta = rand(-pi:pi,2) # Vector de angs aleatorios

        ang = ruido.*eta # aleaorios * eta

        ang_tot = abs(1-w) * (AL[i]+ang[1]) + w * (AS[i]+ang[2])

        # vel_n = RotVec(vel[i],ang_tot)
        vel[i] = RotVec(vel[i],ang_tot)

        # vel[i] = vel_n

    end
end

########################################################

#Actualiza todo

function Evoluciona(i::Int64, step::Int64)

  SR = SetSR(r0,Dist)

  UpdatePos(pos,dt)
  UpdateVel(vel,SR,LR)

  # println(i)

  if  i == 1 || i%step == 0
    println("t = $i writing")
    # println(i)
    PrintTrays(pos)
    PrintDist(i,Dist)
  end

end

########################################################
