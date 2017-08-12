### ============== ### ============== ### ============== ###
##    2D Local and NonLocal Model Simulation Script       ##
##    Martin Zumaya Hernandez                             ##
##    12 / 03 / 2017                                      ##
### ============== ### ============== ### ============== ###

using CollectiveDynamics, Distributions

### ============== ### ============== ###
###          SYSTEM EVOLUTION         ###
### ============== ### ============== ###

function evolve_system(flock, pars, κ_dist, r0)

    ### COMPUTE RELATIVE DISTANCES
    CollectiveDynamics.calc_Rij_MOD(flock.pos, flock.Nij, r0)

    n_nl = rand(κ_dist, pars.N) # obtain non-local input degree for each particle

    ### COMPUTE INTERACTIONS
    # CollectiveDynamics.calc_local_nonLocal_toplogical_interactions!(flock.vel, flock.v_r, flock.v_n, flock.Nij, n_t, n_nl)
    CollectiveDynamics.calc_local_nonLocal_metric_interactions(flock.vel, flock.v_r, flock.v_n, flock.Nij, n_nl)

    ### PARTICLE UPDATE
    map( (p, v, vr, vn) -> CollectiveDynamics.rot_move_part_2D_MOD!(p, v, vr, vn, pars.η, pars.ω), flock.pos, flock.vel, flock.v_r, flock.v_n )

end

### =============== ### =============== ###
###   DEFINITION OF INITIAL PARAMETERS  ###
### =============== ### =============== ###

# N = 128
# κ = 2.0
# ω = 0.5
# T = 3
# rep = 1

N   = parse(Int, ARGS[1]) # number of particles
κ   = parse(Float64, ARGS[2]) # average non-local interactions
ω   = parse(Float64, ARGS[3]) # interactions relative weight
T   = parse(Int, ARGS[4]) # integration time steps
rep = parse(Int, ARGS[5])

ρ = 0.3
η = 0.15

pars = LocNonLocParameters(N, κ, ω, ρ, η) # system's parameters

l = 0.5

L  = sqrt(N / pars.ρ) # size of box

# r0 = (pars.v0 * pars.dt) / pars.l # local interaction range
r0 = (pars.v0 * pars.dt) / l # local interaction range

# p  = pars.κ / (N-1) # non-local link probability

# n_t = 6 # number of topological local interactions
κ_dist = Poisson(κ)

### =============== ### =============== ###
### SET UP SYSTEM AND OUTPUT STRUCTURE  ###
### =============== ### =============== ###

flock = LocNonLocFlock(N, L, pars.v0, 0.0, 2)

output_path = set_output_data_structure_lnl("NLOC_MET_2D", N, κ, ω)

pos_file = open(output_path * "/pos_$(rep).dat", "w+")
vel_file = open(output_path * "/vel_$(rep).dat", "w+")

# net_file = open(output_path * "/net_$(rep).dat", "w+")
# write(net_file, flock.Nij)
# close(net_file)

### ============== ### ============== ###
###          SYSTEM EVOLUTION         ###
### ============== ### ============== ###

# full_time_evolution_2D(pos_file, vel_file, T, flock, r0, pars.η, ω)

times = [convert(Int, exp10(i)) for i in 0:T]

for i in 1:(length(times) - 1)

    if i > 1

        for t in (times[i]+1):times[i+1]

            evolve_system(flock, pars, κ_dist, r0)

            if t % times[i] == 0 || t % times[i-1] == 0
                println("//////// ", t)
                write(pos_file, vcat(flock.pos...))
                write(vel_file, vcat(flock.vel...))
            end
        end

    else

        for t in (times[i]+1):times[i+1]

            evolve_system(flock, pars, κ_dist, r0)

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
