### ============== ### ============== ###
##    3D Simple Vicksek Model          ##
##    Martin Zumaya Hernandez          ##
##    19 / 04 / 2017                   ##
### ============== ### ============== ###

using CollectiveDynamics, Distributions

### ============== ### ============== ###
###          SYSTEM EVOLUTION         ###
### ============== ### ============== ###

function evolve_system(flock, pars, r0, L)

    ### COMPUTE RELATIVE DISTANCES
    CollectiveDynamics.calc_Rij(flock.pos, flock.Nij, r0)

    ### COMPUTE INTERACTIONS
    CollectiveDynamics.calc_vicsek_interactions(flock.vel, flock.v_r, flock.Nij)

    ### PARTICLE UPDATE
    map( (p, v, vr) -> CollectiveDynamics.rot_move_vicsek_3D(p, v, vr, pars.η, L), flock.pos, flock.vel, flock.v_r )

end

### =============== ### =============== ###
###   DEFINITION OF INITIAL PARAMETERS  ###
### =============== ### =============== ###

# N = 128
# κ = 0.0
# ω = 0.5
# T = 3
# rep = 1

N = parse(Int, ARGS[1]) # number of particles
ρ = parse(Float64, ARGS[2]) # density
T = parse(Int, ARGS[3]) # integration time steps
rep = parse(Int, ARGS[4])

η = 0.15

pars = LocNonLocParameters(N, 0.0, 0.0, ρ, η) # system's parameters

L  = cbrt(N / pars.ρ) # size of box

r0 = (pars.v0 * pars.dt) / pars.l # local interaction range
# p  = pars.κ / (N-1) # non-local link probability
# p = 0.0

# n_t = 6 # number of topological local interactions
# κ_dist = Poisson(κ)

# Rij = zeros(Float64, N, N)

### =============== ### =============== ###
### SET UP SYSTEM AND OUTPUT STRUCTURE  ###
### =============== ### =============== ###
flock = LocNonLocFlock(N, L, pars.v0, 0.0, 3)

output_path = CollectiveDynamics.set_output_data_structure_vsk("NLOC_VSK_3D", N, ρ)

pos_file = open(output_path * "/pos_$(rep).dat", "w+")
vel_file = open(output_path * "/vel_$(rep).dat", "w+")

# net_file = open(output_path * "/net_$(rep).dat", "w+")
# write(net_file, flock.Nij)
# close(net_file)

times = [convert(Int, exp10(i)) for i in 0:T]

for i in 1:(length(times) - 1)

    if i > 1

        for t in (times[i]+1):times[i+1]

            evolve_system(flock, pars, r0, L)

            if t % times[i] == 0 || t % times[i-1] == 0
                println("//////// ", t)
                write(pos_file, vcat(flock.pos...))
                write(vel_file, vcat(flock.vel...))
            end
        end

    else

        for t in (times[i]+1):times[i+1]

            evolve_system(flock, pars, r0, L)

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
