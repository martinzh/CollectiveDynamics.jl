module libflock

# using #dependencies

# import #methods to overload

export SetPos, SetVels, SetMatrix, UpdatePos

 #module body


#Genera N vectores aleatorios de dimension dim
function SetPos()
    for i=1:N
        #push!(pos,rand(-L:L,dim)) #inicia en una caja de tamanio L
        push!(pos,rand(dim)*L) #inicia en una caja de tamanio L
        push!(vel,rand(dim))
    end
end

#Normaliza los vectores, magnitud v0 direccion aleatoria
function SetVels()
    for i in 1:N
        vel[i] = vel[i] * (v0/norm(vel[i]))
        #println(norm(vel[i]))
    end
end

#calcula distancias y adjacencia
function SetMatrix()
    # Adj = zeros(N,N) #Limpia matriz adjacencia
    # Dist = zeros(N,N) #Limpia matriz adjacencia
    # Usando sparse
    Adj=spzeros(N,N) #Limpia matriz adjacencia
    Dist=spzeros(N,N) #Limpia matriz adjacencia

    for i=1:N
        for j=N:-1:i
            
            d=norm(pos[i]-pos[j])
            Dist[i;j]=Dist[j;i]=d
            
            #d < r0 && d > 0 ? Adj[i;j] = Adj[j;i] = 1 : Adj[i;j] = 0
            d < r0 ? Adj[i;j]=Adj[j;i]=1 : Adj[i;j]=0
        end
    end
end

#Actualiza posiciones
function UpdatePos()
    #println(pos)
    for i = 1:N
        pos[i] = pos[i] + vel[i] 
    end
#println(pos)
end

end