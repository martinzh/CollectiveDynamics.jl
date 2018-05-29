### ============== ### ============== ### ============== ###
### SHARED ARRAYS VERSION
### METRIC SHORT RANGE INTERACTIONS
### ============== ### ============== ### ============== ###

@everywhere using Distributions, Quaternions
@everywhere include(joinpath(homedir(),"GitRepos","CollectiveDynamics.jl","parallel_test","par_mod.jl"))

### ============== ### ============== ### ============== ###
### SYSTEM'S PARAMETERS
### ============== ### ============== ### ============== ###

file = ARGS[1]

n   = parse(Int64, ARGS[2]) # average non-local interactions
k   = parse(Float64, ARGS[3]) # average non-local interactions
w   = parse(Float64, ARGS[4]) # interactions relative weight
e   = parse(Float64, ARGS[5]) # noise intensity

Ti   = parse(Int, ARGS[6]) # integration time steps
Tf   = parse(Int, ARGS[7]) # integration time steps

rep = parse(Int, ARGS[8])

### ============== ### ============== ### ============== ###

# @everywhere N = 128
# @everywhere κ = 1.0
# @everywhere ω = 0.5
#
# rep = 20

@eval @everywhere N = $n
@eval @everywhere κ = $k
@eval @everywhere ω = $w
@eval @everywhere η = $e

@everywhere ρ = 0.3
@everywhere v0 = 1.0
@everywhere dt = 1.0
@everywhere l = 0.5

@everywhere L  = cbrt(N / ρ) # size of box

@everywhere r0 = ((v0 * dt) / l)^2 # local interaction range

@everywhere κ_dist = Poisson(κ)
# κ_dist = Poisson(κ)

pos = SharedArray{Float64}(3N) # particles positions
vel = SharedArray{Float64}(3N) # array of particles' velocities

v_r = SharedArray{Float64}(3N) # local metric interactions
v_n = SharedArray{Float64}(3N) # non local topological interactions

R_ij = SharedArray{Int64}(N,N)

### ============== ### ============== ### ============== ###

if file == "f"

    ### ============== ### ============== ### ============== ###
    ### INITIALIZATION FROM FILE
    ### ============== ### ============== ### ============== ###

    # pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))
    # vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

    data_path = joinpath(homedir(),"art_DATA","NLOC_MET_3D_EXT_CP","DATA","data_N_$(N)","data_N_$(N)_k_$(κ)_w_$(ω)")

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

    output_path = set_output_data_structure_lnl("NLOC_MET_3D_EXT_LST", N, κ, ω)

    cp(joinpath(data_path,"pos_$(rep).dat"), joinpath(output_path,"pos_$(rep).dat"), remove_destination = true)
    cp(joinpath(data_path,"vel_$(rep).dat"), joinpath(output_path,"vel_$(rep).dat"), remove_destination = true)

    pos_file = open(joinpath(output_path,"pos_$(rep).dat"), "a+")
    vel_file = open(joinpath(output_path,"vel_$(rep).dat"), "a+")

else

    ### ============== ### ============== ### ============== ###
    ### RANDOM INITIAL CONDITIONS
    ### ============== ### ============== ### ============== ###

    for i in 1:length(pos)
        pos[i] = 2*rand()*L - L
        vel[i] = 2*rand() - 1
    end

    for i in 1:3:length(vel)
        norm = sqrt(vel[i]^2 + vel[i+1]^2 + vel[i+2]^2)
        vel[i] /= norm
        vel[i+1] /= norm
        vel[i+2] /= norm
    end

    ### ============== ### ============== ### ============== ###
    ### SET OUTPUT
    ### ============== ### ============== ### ============== ###

    # output_path = set_output_data_structure_lnl("NLOC_MET_3D_EXT", N, κ, ω)
    output_path = set_output_data_structure_lnl("NLOC_MET_3D_EXT", N, κ, ω, η)

    pos_file = open(joinpath(output_path,"pos_$(rep).dat"), "w+")
    vel_file = open(joinpath(output_path,"vel_$(rep).dat"), "w+")

    # write initial conditions
    println("//////// ", 1)
    write(pos_file, pos)
    write(vel_file, vel)
end

### ============== ### ============== ### ============== ###
### TIME EVOLUTION
### ============== ### ============== ### ============== ###

times = [convert(Int, exp10(i)) for i in Ti:Tf]

for i in 1:(length(times) - 1)

    for t in (times[i]+1):times[i+1]

        evolve_metric_system(pos, vel, v_r, v_n, R_ij, r0, N, η, ω, κ_dist)

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
