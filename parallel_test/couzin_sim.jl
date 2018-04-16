### ============ INCLUDE PACKAGES ============ ###

@everywhere using Quaternions
@everywhere include(joinpath(homedir(),"GitRepos","CollectiveDynamics.jl","src","couzin_shared.jl"))

### ============ SYSTEM'S PARAMETERS ============ ###

@everywhere dt = 1.0 # time step
@everywhere ρ = 0.3 # initial density
# @everywhere l = 0.5
@everywhere v0 = 1.0 # speed
@everywhere η = 0.15 # noise intensity
# @everywhere η = 0.25 # noise intensity
@everywhere θ = 40.0 # maximum turn
@everywhere δ = 0.05 # deviation from aligned velocity

### ============ METRIC BEHAVIORAL THRESHOLDS ============ ###

n = parse(Int64, ARGS[1])
# n = 128

o = parse(Float64, ARGS[2])
a = parse(Float64, ARGS[3])
# o = 0.5
# a = 0.5

T = parse(Int64, ARGS[4])

rep = parse(Int64, ARGS[5])

init = ARGS[6] # random or aligned initial velocities

@eval @everywhere init_e = $init

@eval @everywhere N = $n
@eval @everywhere Δo = $o
@eval @everywhere Δa = $a

@everywhere L  = cbrt(N / ρ) # size of box

@everywhere zor = 0.25 # zone of repulsion
# @everywhere zor = 1.0 # zone of repulsion
@everywhere zoo = zor + Δo*L # zone of orientation
@everywhere zoa = zoo + Δa*L # zone of attraction

println("rep:\t", zor, "\torient:\t", zoo, "\tattr:\t", zoa )

### ============ SHARED ARRAY INITIALIZATION ============ ###

pos = SharedArray{Float64}(3N) # particles positions
vel = SharedArray{Float64}(3N) # array of particles' velocities

v_int = SharedArray{Float64}(3N) # total signal
v_r   = SharedArray{Float64}(3N) # repulsion metric interactions
v_o   = SharedArray{Float64}(3N) # orientation metric interactions
v_a   = SharedArray{Float64}(3N) # attraction metric interactions

Rij = SharedArray{Float64}(N,N)

output_path = ""

if init_e == "R"

    output_path = set_output_data_structure_lnl("COUZIN_3D_TEST", N, ARGS[2], ARGS[3])
    println(output_path)

    ### ============ RANDOM INITIAL CONDITIONS ============ ###
    for i in 1:length(pos)
        pos[i] = 2*rand()*L - L
        vel[i] = 2*rand() - 1
    end

elseif init_e == "A"

    output_path = set_output_data_structure_lnl("COUZIN_3D_VAL", N, ARGS[2], ARGS[3])
    println(output_path)

    ### ============ RANDOM POSITIONS BUT ALIGNED VELOCITIES ============ ###

    vel_0 = [2*rand() - 1, 2*rand() - 1, 2*rand() - 1]

    for i in 1:3:length(pos)
        pos[i]   = 2*rand()*L - L
        pos[i+1] = 2*rand()*L - L
        pos[i+2] = 2*rand()*L - L
        vel[i]   = vel_0[1] + δ*rand() - δ
        vel[i+1] = vel_0[2] + δ*rand() - δ
        vel[i+2] = vel_0[3] + δ*rand() - δ
    end

end

### ============ VELOCITY NORMALIZATION ============ ###

for i in 1:3:length(vel)
    norm = sqrt(vel[i]^2 + vel[i+1]^2 + vel[i+2]^2)
    vel[i] /= norm
    vel[i+1] /= norm
    vel[i+2] /= norm
end

### ============ SET OUTPUT ============ ###

pos_file = open(joinpath(output_path,"pos_$(rep).dat"), "w+")
vel_file = open(joinpath(output_path,"vel_$(rep).dat"), "w+")

# write initial conditions
write(pos_file, pos)
write(vel_file, vel)
println("//////// ", 1)

### ============ TIME EVOLUTION ============ ###

times = [convert(Int, exp10(i)) for i in 0:T]

for i in 1:(length(times) - 1)

    for t in (times[i]+1):times[i+1]

        evolve_system(pos, vel, v_int, v_r, v_o, v_a, Rij, N, zor, zoo, zoa, η, θ, v0)
        # evolve_system_vr(pos, vel, v_r, Rij, N, zor, zoo, zoa, η, θ, v0)

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
