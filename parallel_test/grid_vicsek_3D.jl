### ============== ### ============== ###
##    3D Simple Vicksek Model          ##
##    Grid version                     ##
##    (cell velocity search)           ##
##    Martin Zumaya Hernandez          ##
##    30 / 06 / 2017                   ##
### ============== ### ============== ###

### ============== ### ============== ###

using CollectiveDynamics

include("$(homedir())/GitRepos/CollectiveDynamics.jl/parallel_test/SVM.jl")

### ============== ### ============== ###

### ============== ### ============== ###
###          SYSTEM EVOLUTION         ###
### ============== ### ============== ###

function evolve_system(flock, box, cell_size)

    # compute each particle cell_id and assgin particles id to each box
    assign_cell(flock, box, cell_size)

    # for i in eachindex(flock.v_r)
    #     println(i, "\t",flock.vel[i])
    # end

    # Search interactions in adjacent cells
    for i in eachindex(flock.pos)
        get_neighbors(i, flock, box)
        # println("pass")
    end

    # for i in eachindex(flock.v_r)
    #     println(i, "\t",flock.v_r[i])
    # end

    # update particles' position
    broadcast(update_part, flock.pos, flock.vel, flock.v_r, flock.η, box.L)

end

### =============== ### =============== ###
###   DEFINITION OF INITIAL PARAMETERS  ###
### =============== ### =============== ###

# ρ = parse(Float64, ARGS[1]) # densityz
# T = parse(Int, ARGS[2]) # integration time steps
# rep = parse(Int, ARGS[3])

### =============== ### =============== ###

N = parse(Int, ARGS[1]) # number of particles
ρ = parse(Float64, ARGS[2]) # density
T = parse(Int, ARGS[3]) # integration time steps
rep = parse(Int, ARGS[4])

# N = 2^13
# N = 2^11
N = 4096

ρ = 0.1
# T = 3
# rep = 1

η = 0.15

v0 = 1.0
dt = 1.0

### =============== ### =============== ###

M = 10 # number of boxes per dimension
# M = 5 # number of boxes per dimension
# M = 4 # number of boxes per dimension

l = 0.5 # interaction range is double the distance each particle moves in one time step
# l = 0.1 # interaction range is ten times the distance each particle moves in one time step
r0 = (v0 * dt) / l # local interaction range

# L = M * r0 # box size in terms of interaction range
L  = cbrt(N / ρ) # size of box

### =============== ### =============== ###

# N = convert(Int, ρ * (M * r0)^3)
# println("N = ",N)

cell_size = step(linspace(0., L , M))

### ============== ### ============== ###

box = Box(L, M)
flock = Flock(N, L, dt, v0, r0, η)

### ============== ### ============== ###

output_path = CollectiveDynamics.set_output_data_structure_vsk("SVM_GRID_FN_3D", N, ρ)

pos_file = open(output_path * "/pos_$(rep).dat", "w+")
vel_file = open(output_path * "/vel_$(rep).dat", "w+")

times = [convert(Int, exp10(i)) for i in 0:T]

for i in 1:(length(times) - 1)

    if i > 1

        for t in (times[i]+1):times[i+1]

            evolve_system(flock, box, cell_size)

            if t % times[i] == 0 || t % times[i-1] == 0
                println("//////// ", t)
                write(pos_file, vcat(flock.pos...))
                write(vel_file, vcat(flock.vel...))
            end
        end

    else

        for t in (times[i]+1):times[i+1]

            evolve_system(flock, box, cell_size)

            if t % times[i] == 0
                println("//////// ", t)
                write(pos_file, vcat(flock.pos...))
                write(vel_file, vcat(flock.vel...))
            end
        end

    end

end

close(pos_file)
close(vel_file)

println("Done all")


# for i in 1:1000
#
#     println(i)
#
#     assign_cell(flock, box)
#
#     # for i in eachindex(flock.v_r)
#     #     println(i, "\t",flock.vel[i])
#     # end
#
#     # Search interactions in adjacent cells
#     for i in eachindex(flock.pos)
#         get_neighbors(i, flock, box)
#         # println("pass")
#     end
#
#     # for i in eachindex(flock.v_r)
#     #     println(i, "\t",flock.v_r[i])
#     # end
#
#     # update particles' position
#     broadcast(update_part, flock.pos, flock.vel, flock.v_r, flock.η, box.L)
#
# end
