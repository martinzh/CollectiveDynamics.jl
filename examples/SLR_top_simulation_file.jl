### ============== ### ============== ###
##    Short and Long-range             ##
##    interactions models              ##
##    (topological short-range)        ##
##    Martin Zumaya Hernandez          ##
##    EXAMPLE SIMULATION SCRIPT        ##
### ============== ### ============== ###

@everywhere using Distributions, CollectiveDynamics.ShortLongRange

### ============== ### ============== ### ============== ###
### SYSTEM'S PARAMETERS
### ============== ### ============== ### ============== ###

n   = parse(Int64, ARGS[1]) # number of particles
k   = parse(Float64, ARGS[2]) # average non-local interactions
w   = parse(Float64, ARGS[3]) # interactions relative weight
e   = parse(Float64, ARGS[4]) # noise intensity

Ti   = parse(Int, ARGS[5]) # start of integration time steps (10^Ti)
Tf   = parse(Int, ARGS[6]) # end of integration time steps (10^Tf)

rep = parse(Int, ARGS[7])

### ============== ### ============== ### ============== ###

@eval @everywhere N = $n
@eval @everywhere κ = $k
@eval @everywhere ω = $w
@eval @everywhere η = $e

@everywhere k_sh = 6

@everywhere ρ = 0.3
@everywhere v0 = 1.0
@everywhere dt = 1.0
@everywhere l = 0.5

@everywhere L  = cbrt(N / ρ) # size of box

@everywhere κ_dist = Poisson(κ)

pos = SharedArray{Float64}(3N) # particles positions
vel = SharedArray{Float64}(3N) # array of particles' velocities

v_r = SharedArray{Float64}(3N) # local metric interactions
v_n = SharedArray{Float64}(3N) # non local topological interactions

R_ij = SharedArray{Float64}(N,N)

### ============== ### ============== ### ============== ###
### INITIALIZATION FROM FILE
### ============== ### ============== ### ============== ###

# pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))
# vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

data_path = joinpath(homedir(),"art_DATA","testing", "NLOC_TOP_3D_EXT_CP","DATA","data_N_$(N)","data_N_$(N)_k_$(κ)_w_$(ω)")

raw_data = reinterpret(Float64, read(joinpath(data_path,"pos_$(rep).dat")))
pos_data = raw_data[(end-3N+1):end]

raw_data = reinterpret(Float64, read(joinpath(data_path,"vel_$(rep).dat")))
vel_data = raw_data[(end-3N+1):end]

for i in 1:length(pos)
    pos[i] = pos_data[i]
    vel[i] = vel_data[i]
end


### ============== ### ============== ### ============== ###
### SET OUTPUT
### ============== ### ============== ### ============== ###

output_path = set_output_data_structure("SLR_TOP", N, κ, ω)
# output_path = set_output_data_structure("SLR_TOP", N, κ, ω, η)

pos_file = open(joinpath(output_path,"pos_$(rep).dat"), "w+")
vel_file = open(joinpath(output_path,"vel_$(rep).dat"), "w+")

# write initial conditions
println("//////// ", 1)
write(pos_file, pos)
write(vel_file, vel)

### ============== ### ============== ### ============== ###
### TIME EVOLUTION
### ============== ### ============== ### ============== ###

times = [convert(Int, exp10(i)) for i in Ti:Tf]

for i in 1:(length(times) - 1)

    for t in (times[i]+1):times[i+1]

        evolve_topological_system(pos, vel, v_r, v_n, R_ij, N, k_sh, η, ω, κ_dist)

        if t % times[i] == 0 || t % div(times[i], exp10(1)) == 0
            println("//////// ", t)
            write(pos_file, pos)
            write(vel_file, vel)
        end
    end

end

close(pos_file)
close(vel_file)

rmprocs(workers())

println("Done all")

### ============== ### ============== ### ============== ###
