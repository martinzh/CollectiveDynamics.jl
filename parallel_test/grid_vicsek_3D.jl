### ============== ### ============== ###
##    3D Simple Vicksek Model          ##
##    Grid version                     ##
##    Martin Zumaya Hernandez          ##
##    30 / 06 / 2017                   ##
### ============== ### ============== ###

# using CollectiveDynamics, Distributions
using Plots, CollectiveDynamics
using CollectiveDynamics
gr()

### ============== ### ============== ###
###          SYSTEM EVOLUTION         ###
### ============== ### ============== ###

### =============== ### =============== ###
###   DEFINITION OF INITIAL PARAMETERS  ###
### =============== ### =============== ###

# N = parse(Int, ARGS[1]) # number of particles
# ρ = parse(Float64, ARGS[2]) # density
# T = parse(Int, ARGS[3]) # integration time steps
# rep = parse(Int, ARGS[4])

N = 2^13
N = 2^10


T = 3
rep = 1

η = 0.15

pars = LocNonLocParameters(N, 0.0, 0.0, ρ, η) # system's parameters

ρ = 0.3

L  = cbrt(N / ρ) # size of box

l = 1.0
r0 = (pars.v0 * pars.dt) / pars.l # local interaction range
r0 = (pars.v0 * pars.dt) / l # local interaction range

flock = LocNonLocFlock(N, L, pars.v0, 0.0, 3)

M = 10. # number of boxes per dimension

b_id = [zeros(Int,3) for i in 1:N] # box index for each particle in each dimension

pos = [ [L*rand(), L*rand(), L*rand()] for i in 1:N ]

### ============== ### ============== ###

x = [pos[i][1] for i in 1:N]
y = [pos[i][2] for i in 1:N]
z = [pos[i][3] for i in 1:N]

scatter(x,y,z, alpha = 0.4, leg = false, ms = 3)

### ============== ### ============== ###

for i in 1:3
    println(pos[1][i],"\t",div(floor.(pos[1][i]),L/M))
end

p_num = Dict()

for i in 1:N
    # b_id[i] = convert(Array{Int}, div.(floor.(pos[i]), L/M)) + 1
    b_id[i] = convert(Array{Int}, div.(floor.(pos[i]), r0)) + 1
    haskey(p_num, b_id[i]) ? p_num[b_id[i]] += 1 : p_num[b_id[i]] = 1
end

x = [pos[i][1] for i in find(x -> x==[1,1,1], b_id)]
y = [pos[i][2] for i in find(x -> x==[1,1,1], b_id)]
z = [pos[i][3] for i in find(x -> x==[1,1,1], b_id)]

convert(Array{Int}, div.(floor.(pos[1]), L/M) + 1)

maximum(values(p_num))

findmax(values(p_num))
