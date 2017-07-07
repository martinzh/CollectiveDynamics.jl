### ============== ### ============== ###
##    3D Simple Vicksek Model          ##
##    Grid version (particle search)   ##
##    Martin Zumaya Hernandez          ##
##    30 / 06 / 2017                   ##
### ============== ### ============== ###

# using CollectiveDynamics, Distributions
using Plots
using CollectiveDynamics
gr()
pyplot()

gui()

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

# N = 2^13
N = 2^10

ρ = 0.5
T = 3
rep = 1

η = 0.15

v0 = 1.0
dt = 1.0

### =============== ### =============== ###

M = 10 # number of boxes per dimension
M = 5 # number of boxes per dimension
M = 4 # number of boxes per dimension

l = 0.5 # interaction range is double the distance each particle moves in one time step
l = 0.1 # interaction range is ten times the distance each particle moves in one time step
r0 = (v0 * dt) / l # local interaction range

L = M * r0 # box size in terms of interaction range

L  = cbrt(N / ρ) # size of box

### =============== ### =============== ###

N = convert(Int, ρ * (M * r0)^3)

### ============== ### ============== ###

bulk_cells = [[i,j,k] for i in 2:M-1, j in 2:M-1, k in 2:M-1]

bottom_cells = union([[i,j,1] for i in 2:M-1, j in 1:M], [[i,j,1] for i in [1,M], j in 2:M-1])
top_cells = union([[i,j,M] for i in 2:M-1, j in 1:M], [[i,j,1] for i in [1,M], j in 2:M-1])

left_wall = union([[1,j,k] for j in 2:M-1, k in 1:M], [[1,j,k] for j in [1,M], k in 2:M-1])
right_wall = union([[M,j,k] for j in 2:M-1, k in 1:M], [[M,j,k] for j in [1,M], k in 2:M-1])

front_wall = union([[i,1,k] for i in 2:M-1, k in 1:M], [[i,1,k] for i in [1,M], k in 2:M-1])
back_wall = union([[i,M,k] for i in 2:M-1, k in 1:M], [[i,M,k] for i in [1,M], k in 2:M-1])

corners = [[i,j,k] for i in [1,M], j in [1, M], k in [1,M]]
### ============== ### ============== ###

b_id = [zeros(Int,3) for i in 1:N] # box index for each particle in each dimension

pos = [ [L*rand(), L*rand(), L*rand()] for i in 1:N ]
vel = v0 * [ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:N ]

part_box = Dict{Array{Int,1}, Array{Int, 1}}() # particles id's in each box

vel_box = Array{Array{Float64,1}}(M,M,M) # average velocity per box
center_box = Array{Array{Float64,1}}(M,M,M) # position of center of each box
# center_box = Array{Array{Float64,1}}(M-1,M-1,M-1) # position of center of each box

### ============== ### ============== ###
# Initialize velocity vector perop box and compute center of each box

for i in 1:length(vel_box)
    vel_box[i] = zeros(3)
end

r = linspace(0., L, M+1)
mp = collect(r[1:length(r) - 1] + 0.5 * step(r))

for i in 1:M, j in 1:M, k in 1:M
    # println(i,"\t", j, "\t", k)
    center_box[i,j,k] = [mp[i], mp[j], mp[k]]
end

### ============== ### ============== ###

# compute each particle box_id and assgin particles id to each box

for i in 1:N
    b_id[i] = convert(Array{Int}, div.(floor.(pos[i]), L/M)) + 1
    # b_id[i] = convert(Array{Int}, div.(floor.(pos[i]), r0)) + 1

    haskey(part_box, b_id[i]) ? push!(part_box[b_id[i]], i) : part_box[b_id[i]] = [i]
end

for kv in part_box
    println(kv)
end
### ============== ### ============== ###

# compute each box average velocity

for k in keys(part_box)
    vel_box[k[1], k[2], k[3]] = mean([vel[i] for i in part_box[k]])
end

### ============== ### ============== ###
## Search interactions in adjacent cells
i = 1
b_id[i]


### ============== ### ============== ###

x = [pos[i][1] for i in 1:N]
y = [pos[i][2] for i in 1:N]
z = [pos[i][3] for i in 1:N]

scatter(x,y,z, alpha = 0.4, leg = false, ms = 3)

x = [pos[i][1] for i in find(x -> x==[1,1,1], b_id)]
y = [pos[i][2] for i in find(x -> x==[1,1,1], b_id)]
z = [pos[i][3] for i in find(x -> x==[1,1,1], b_id)]

### ============== ### ============== ###


img = plot()

for j in 1:M-1

    x = hcat(center_box...)[1,:]
    y = hcat(center_box...)[2,:]
    z = hcat(center_box...)[3,:]

    scatter!(x, y, z, alpha = 0.8, leg = false, ms = 2)

end


### ============== ### ============== ###

using PyPlot

pygui(true)

fig = figure()
ax = fig[:gca](projection="3d")

x = hcat(center_box...)[1,:]
y = hcat(center_box...)[2,:]
z = hcat(center_box...)[3,:]

u = 2*hcat(vel_box...)[1,:]
v = 2*hcat(vel_box...)[2,:]
w = 2*hcat(vel_box...)[3,:]

ax[:quiver](x,y,z, u,v,w)

### ============== ### ============== ###
# search for interacting particles in bulk cells

c = 0

println("i\tj\tk\tc")

for i in -1:1, j in -1:1, k in -1:1
    c += 1
    println(i, "\t", j, "\t", k, "\t", c)
end

### ============== ### ============== ###
# search for interacting particles in boundary cells

c = 0

println("i\tj\tk\tc")

for i in -1:1, j in -1:1, k in -1:1
    c += 1
    println(i, "\t", j, "\t", k, "\t", c)
end

### ============== ### ============== ###
# bulk_cells

in(3, 2:4)

in(b_id[1], bulk_cells)
