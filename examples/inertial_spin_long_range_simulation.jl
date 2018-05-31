### ============== ### ============ ###
##    Extended Inertial Spin Model   ##
##    Martin Zumaya Hernandez        ##
##    EXAMPLE SIMULATION SCRIPT      ##
### ============== ### ============ ###

### ============ INCLUDE PACKAGES ============ ###

using CollectiveDynamics.InertialSpin

### =============== ### =============== ###
###   DEFINITION OF INITIAL PARAMETERS  ###
### =============== ### =============== ###

# N   = 128
# η   = 1.0
# τ   = 6
# rep = 1
# T   = 8*exp10(-5)
# δ   = 0.35

N    = parse(Int, ARGS[1]) # Number of particles
η    = parse(Float64, ARGS[2]) # Dissipation term
T    = parse(Float64, ARGS[3]) # temperature, noise
n_nl = parse(Float64, ARGS[4]) # non local interactions per particle
τ    = parse(Int, ARGS[5]) # number of iterations
rep  = parse(Int, ARGS[6]) # repetition number

χ   = 1.25
J   = 0.8
ρ   = 0.3
v0  = 0.1
n_t = 6

pars = InertialParameters(χ, J, η, v0*sqrt(J/χ), ρ, v0, N, T, n_t, n_nl)

σ = sqrt((2pars.d) * η * T) # noise std deviation ( square root of variance )
L = cbrt(N / pars.ρ) # 3D

### =============== ### =============== ###
### SET UP SYSTEM AND OUTPUT STRUCTURE  ###
### =============== ### =============== ###

flock = InertialExtFlock(N, L, pars.v0)

output_path = set_output_data_structure_lr("EXTENDED_INERTIAL_SPIN", N, η, T, n_nl)

pos_file  = open(output_path * "/pos_$(rep).dat", "w+")
vel_file  = open(output_path * "/vel_$(rep).dat", "w+")
spin_file = open(output_path * "/spin_$(rep).dat", "w+")

### ============== ### ============== ###
###          SYSTEM EVOLUTION         ###
### ============== ### ============== ###

times = [convert(Int, exp10(i)) for i in 0:τ]

  for i in 1:(length(times) - 1)

      if i > 1

          for t in (times[i]+1):times[i+1]

              evolve_extended_system(flock.pos, flock.vel, flock.v_t, flock.v_nl, flock.spin, flock.Rij, pars, σ)

              if t % times[i] == 0 || t % times[i-1] == 0
                  println("//////// ", t)
                  write(pos_file, vcat(flock.pos...))
                  write(vel_file, vcat(flock.vel...))
                  write(spin_file, vcat(flock.spin...))
              end
          end

      else

          for t in (times[i]+1):times[i+1]

              evolve_extended_system(flock.pos, flock.vel, flock.v_t, flock.v_nl, flock.spin, flock.Rij, pars, σ)

              if t % times[i] == 0
                  println("//////// ", t)
                  write(pos_file, vcat(flock.pos...))
                  write(vel_file, vcat(flock.vel...))
                  write(spin_file, vcat(flock.spin...))
              end
          end

      end

  end

close(pos_file)
close(vel_file)

println("Done all")
