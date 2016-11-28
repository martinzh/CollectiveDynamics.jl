### ============== ### ============== ### ============== ###
## Numerical simulations of the Discretized version of the
## inertial spin model Reference: Cavagna...
## Martin Zumaya Hernandez
## 14 / 11 / 2016
### ============== ### ============== ### ============== ###

using CollectiveDynamics

### ============== ### ============== ### ============== ###
### PARAMETERS ##
### ============== ### ============== ### ============== ###

immutable parameters
    N  ::Real  # Number of particles
    κ  ::Real # Connectivity
    η  ::Real  # generalized temperature
    ρ  ::Real # Local initial density
    ω  ::Real # Interactions relative weight
    v0 ::Real  # Speed
    r0 ::Real  # Local interaction range
    dt ::Real  # Integration step

    # default constructor
    parameters() = new(1024, 0.0, 0.15, 0.3, 1.0, 10.0, 1.0)

    # full constructor
    parameters(N, κ, η, ρ, ω, v0, r0) = new(N, κ, η, ρ, ω, v0, r0, 1.0)
end

### ============== ### ============== ### ============== ###
### SPIN UPDATE
### ============== ### ============== ### ============== ###

function local_int(pos, vel, v_r, rij, v_n, pars)

    # # matrix of relative distances between particles
    # nij = zeros(pars.N,pars.N)

    # compute nij entries
    for i in 1:pars.N, j = i:pars.N
        rij[i,j] = norm(pos[i] - pos[j])
        rij[j,i] = rij[i,j]
    end

    # compute local metric interaction
    for i in 1:pars.N
        v_r[i] = mean([vel[i] for i in find(x -> x <= pars.r0 && x > 0.0, rij[:,i])])
    end

    # compute non local topological interaction
    for i in 1:pars.N
        v_n[i] = mean([vel[i] for i in find(x -> x <= pars.r0 && x > 0.0, rij[:,i])])
    end



end

### ============== ### ============== ### ============== ###
### SYSTEM EVOLUTION
### ============== ### ============== ### ============== ###

function evolve(pos, vel, vn, spin, pars, σ)

    ### VELOCITY UPDATE
    map!((x,y) -> y + (1/pars.χ)*cross(x,y), vel, spin, vel)

    ### SPIN UPDATE
    spin_update(pos, vel, v_n, spin, pars)

    ### POSITION UPDATE
    map!( (x,y) -> x + y*pars.dt, pos, pos, vel )

end


### =============== ### =============== ### =============== ###
### DEFINITION OF INITIAL PARAMETERS
### =============== ### =============== ### =============== ###

N  = 512 # number of particles
κ = 7 # average non-local interactions
η = 0.15 # noise intensity
ρ  = 0.3 # density
L  = sqrt(N / ρ) # size of box

l = 0.1 # Regimen de velocidad
dt = 1.0 # Time integration step

v0 = 1.0 # particle's speed
T  = 100 # integration time steps

r0 = (v0 * dt) / l

p = κ / (N-1)

pars = parameters( N, κ, η, ρ, 0.5, v0, r0 )

### =============== ### =============== ### =============== ###
### INITIALIZATION OF PARTICLES INITIAL POSITIONS AND VELOCIDITES
### =============== ### =============== ### =============== ###

# array of random initial particles' postitions
# pos = [ [2*rand()*L - L, 2*rand()*L - L, 2*rand()*L - L] for i in 1:pars.N ]
pos = [ [2*rand()*L - L, 2*rand()*L - L] for i in 1:pars.N ]

# array of particles' velocities
# vel = pars.v0 * [ normalize([2*rand() - 1, 2*rand() - 1, 2*rand() - 1]) for i in 1:pars.N ]
vel = pars.v0 * [ normalize([2*rand() - 1, 2*rand() - 1]) for i in 1:pars.N ]

# local metric interactions
# v_n = [zeros(3) for i in 1:pars.N]
v_r = [zeros(2) for i in 1:pars.N]

# non local topological interactions
# v_n = [zeros(3) for i in 1:pars.N]
v_n = [zeros(2) for i in 1:pars.N]

rij = zeros(pars.N, pars.N)

Nij = zeros(pars.N, pars.N)

for i in 1:N, j in union(1:i-1, i+1:N)
    rand() < κ / (N-1) ? Nij[j, i] = 1 : Nij[j, i] = 0
end

### ============== ### ============== ### ============== ###
### SYSTEM EVOLUTION
### ============== ### ============== ### ============== ###

for t in 1:τ
    println(t, "\t", pos[1])
    evolve(pos, vel, v_n, spin, pars, σ)
end
