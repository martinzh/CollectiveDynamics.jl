### ============== ### ============== ###
##    3D Simple Vicksek Model          ##
##    Grid version (particle search)   ##
##    Martin Zumaya Hernandez          ##
##    30 / 06 / 2017                   ##
### ============== ### ============== ###

### ============== ### ============== ###

using CollectiveDynamics.SVM3D

### =============== ### =============== ###
###   DEFINITION OF INITIAL PARAMETERS  ###
### =============== ### =============== ###

N   = parse(Int, ARGS[1]) # number of particles
ρ   = parse(Float64, ARGS[2]) # density
η   = parse(Float64, ARGS[3]) # noise intensity
T   = parse(Int, ARGS[4]) # integration time steps
rep = parse(Int, ARGS[5]) # ensemble index

### =============== ### =============== ###

v0 = 1.0 # particle's speed
dt = 1.0 # integration time step

### =============== ### =============== ###

M = 10 # number of boxes per dimension

l = 0.5 # interaction range is double the distance each particle moves in one time step
# l = 0.1 # interaction range is ten times the distance each particle moves in one time step

r0 = (v0 * dt) / l # local interaction range

L  = cbrt(N / ρ) # size of box

### =============== ### =============== ###

cell_size = step(linspace(0., L , M))

### ============== ### ============== ###

box = Box(L, M)
flock = Flock(N, L, dt, v0, r0, η)

### ============== ### ============== ###

output_path = set_output_data_structure_vsk("SVM_3D", N, ρ)

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
