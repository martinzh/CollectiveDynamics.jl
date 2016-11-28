### ============== ### ============== ### ============== ###
## Numerical simulations of the Discretized version of the
## inertial spin model Reference: Cavagna...
## Martin Zumaya Hernandez
## 14 / 11 / 2016
### ============== ### ============== ### ============== ###

# using CollectiveDynamics

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

function calc_interactions(pos, vel, v_r, rij, v_n, Nij, pars)

    # compute nij entries
    for i in 1:pars.N, j = i:pars.N
        rij[i,j] = norm(pos[i] - pos[j])
        rij[j,i] = rij[i,j]
    end

    # compute local and non local interaction
    for i in 1:pars.N
        v_r[i] = mean([vel[i] for i in find(x -> x <= pars.r0 && x > 0.0, rij[:,i])])

        if length(findn(Nij[:,i])) != 0
            v_n[i] = mean([vel[i] for i in findn(Nij[:,i])])
        else
            v_n[i] = zeros(2)
        end

    end
end

function rot_move(pos, vel, v_r, v_n, pars)

    for i in 1:pars.N

        prop_angle = atan2(vel[i][2], vel[i][1])

        i_vx = v_r[i][1]
        i_vy = v_r[i][2]

        if i_vx != 0.0 || i_vy != 0.0
            loc_angle = atan2(i_vy, i_vx) - prop_angle;
        else
            loc_angle = 0.0;
        end

        i_vx = v_n[i][1]
        i_vy = v_n[i][2]

        if i_vx != 0.0 || i_vy != 0.0
            non_loc_angle = atan2(i_vy, i_vx) - prop_angle;
        else
            non_loc_angle = 0.0;
        end

        total_angle = pars.ω * loc_angle + (1 - pars.ω) * non_loc_angle + pars.η * (2*rand()*pi - pi);

        c = cos(total_angle)
        s = sin(total_angle)

        vx = vel[i][1]*c - vel[i][2]*s;
        vy = vel[i][1]*s + vel[i][2]*c;

        vel[i][1] = vx;
        vel[i][2] = vy;

        pos[i][1] += vx;
        pos[i][2] += vy;

    end

end

### ============== ### ============== ### ============== ###
### SYSTEM EVOLUTION
### ============== ### ============== ### ============== ###

function evolve(pos, vel, v_r, v_n, rij, Nij, pars)

    calc_interactions(pos, vel, v_r, rij, v_n, Nij, pars)
    rot_move(pos, vel, v_r, v_n, pars)

end


### =============== ### =============== ### =============== ###
### DEFINITION OF INITIAL PARAMETERS
### =============== ### =============== ### =============== ###

N  = 1024 # number of particles
κ = 7 # average non-local interactions
η = 0.15 # noise intensity
ρ  = 0.3 # density
L  = sqrt(N / ρ) # size of box

l = 0.1 # Regimen de velocidad
dt = 1.0 # Time integration step

v0 = 1.0 # particle's speed
T  = 10 # integration time steps

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

for t in 1:T
    println(t, "\t", pos[1])
    evolve(pos, vel, v_r, v_n, rij, Nij, pars)
end
