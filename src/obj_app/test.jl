include("obj_lib.jl")
using DataFrames
## =========================== ## ## =========================== ##


T = 10
step = 1
dt = 1.0
v0 = 1.0
ρ = 10.0
L = 1.0
l = 0.25

N = int(L*L*ρ)

r0 = v0 * dt / l

k = 4
## =========================== ## ## =========================== ##

# Si ω = 0 → solo interaciiones de la red
# Si ω = 1 → solo interaciiones geometricas

ω = 0.0

## =========================== ## ## =========================== ##

# Intensidad de ruido

η = 0.0

## =========================== ## ## =========================== ##
#RandVec(L)
#RandNum(L)

# for i in 1:5
#   println("$(Inputs(k,N,1)),$i")
# end

# Inputs(k,N,2)
## =========================== ## ## =========================== ##
parts = Array(Bird,N)
InitParts(N,L,v0,k)

parts[1].vel

v = parts[1].vel
v += [1.0,1.0]
v
for i in parts[1].inputs
  println(i)
end

GetAngsIN(parts)

Dist = zeros(N,N)

SetSR(r0,Dist,parts)

UpdatePos!(parts,dt)


GetAngs(parts,SetSR(r0,Dist,parts))
GetAngsIN(parts)
