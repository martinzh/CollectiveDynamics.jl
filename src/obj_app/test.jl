include("obj_lib.jl")

N  = 5
L  = 10.0
v0 = 5.0
dt = 1.0
k = 2
r0 = 5.0
w = 0.15
η = 0.1
# una = Bird(RandVec(1.0),RandVec(1.0))
# una = Bird(RandVec(5.0) , RandVec(5.0))

# println(una.pos)
# println(typeof(una))

parts = Array(Bird,N)

typeof(parts)
# InitParts()
println("Inicializa:")
@time InitParts()

for bird = parts
  println(bird.pos)
end

typeof(parts[1].vel)

vel = RandVec(1.0)

typeof(vel)

scale!(vel,v0/norm(vel))

typeof(vel)

norm(vel)

s = [1.0,2.0,0.0]

for i = s
  println(i)
end

typeof(s)

println(size(parts,1))

Dist = zeros(N,N) #Matriz de distancias

#Usando sparse
LR = spzeros(N,N) #Interacciones de largo alcanze
                  #No cambia en el tiempo

# SetLR(k,LR)
println("Crea no-contacto:")
@time SetLR(k,LR)

LR
# println(LR)
# println(issparse(LR))

# rows = rowvals(LR)
# # vals = nonzeros(LR)
# m, n = size(LR)

# println(rows)
# # println(vals)

# for i = 1:n

# 	# println(nzrange(LR,i))
# 	s = 0

#    for j in nzrange(LR, i)
#    		k = rows[j]
#    		println("$i\t$k")
#    		s += 1
#    #    row = rows[j]
#    #    val = vals[j]

#    #    println("$row\t$j")

#    #    # perform sparse wizardry...
#    end
#    println("k = $s")
# end

# println(GetAngs(parts,LR))
println("GetAngs:")
@time GetAngs(parts,LR)

# for i = 1:size(LR)[1]
# 	println(i)
# 	# println(LR[:;i])
# 	# println(LR[1:end,i])

# 	println()
# end


println("Calcula dist y adj:")
@time SR = SetSR(r0,Dist,parts)

# println("Distancias:")
# println(Dist)

# println("Adjacencias:")
# println(SR)

# println(LR)

# println(size(parts)[1])

# println(typeof(parts))
# println(parts)

# for i = 1:N
# 	norma = norm(parts[i].vel)
# 	println(norma)
# end

# for bird = parts
# 	println(bird.pos)
# end

# println("update")

# UpdatePos(parts,dt)
println("Actualiza posicion:")
@time UpdatePos!(parts,dt)

println("Actualiza velocidad:")
@time UpdateVel!(parts,SR,LR,η,w)

# for bird = parts
# 	println(repr(bird.pos)[2:end-1])
# end

V = [1.0 0.0]
α = pi/2

@time RotVec(V,α)
@time RotVec!(V,α)
V
size(parts,1)
RotVec!(parts[1].vel,α)

@time RandVec(5.0)
@time RandVec1(5.0)
