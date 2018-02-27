
using Quaternions

### ============ SYSTEM'S PARAMETERS ============ ###

N = 10

ρ = 0.3
η = 0.15
v0 = 1.0
dt = 1.0
l = 0.5

L  = cbrt(N / ρ) # size of box

### ============ METRIC BEHAVIORAL THRESHOLDS ============ ###

zor = ((v0 * dt) / l)^2 # zone of repulsion
zoo = zor + ((v0 * dt) / l)^2 # zone of orientation
zoa = zoo + ((v0 * dt) / l)^2 # zone of attraction

### ============ RANDOM INITIAL CONDITIONS ============ ###

pos = [[-L + 2*rand()*L, -L + 2*rand()*L, -L + 2*rand()*L] for i in 1:N] # particles positions
vel = [[-1.0 + 2*rand(), -1.0 + 2*rand(), -1.0 + 2*rand()] for i in 1:N] # array of particles' velocities

v_r = [zeros(Float64, 3) for i in 1:N]

for i in 1:N normalize!(vel[i]) end

### ============ COMPUTE RELATIVE DISTANCES ============ ###

Rij = zeros(Float64, N, N)

for j in 1:N, i in 1+j:N
    Rij[i,j] = norm(pos[i] - pos[j])
end

Rij = Symmetric(Rij, :L)

### ============ COMPUTE INTERACTIONS ============ ###

for i in 1:N

    println(i)

    v_r[i] = zeros(Float64, 3)

    repel_neighbors = find( x-> x > 0.0 && x <= zor, Rij[:, i])

    if length(repel_neighbors) > 0

        for j in repel_neighbors
            v_r[i] -= pos[j] - pos[i] / Rij[j,i]
        end

    else

        orient_neighbors = find( x-> x > zor && x < zoo, Rij[:, i])
        atract_neighbors = find( x-> x > zoo && x < zoa, Rij[:, i])

        v_o = zeros(Float64, 3)
        v_a = zeros(Float64, 3)

        for j in find( x-> x > zor && x < zoo, Rij[:, i])
            v_o += vel[j]
        end

        for j in find( x-> x > zoo && x < zoa, Rij[:, i])
            v_a += pos[j] - pos[i] / Rij[j,i]
        end

        if length(v_a) > 0
            v_r[i] = 0.5 * (v_o + v_a)
        else
            v_r[i] = v_o
        end

    end

end

### ============ UPDATE PARTICLE'S POSITIONS AND VELOCITIES ============ ###

for i in 1:N

    println(i)

    q_r = Quaternion(zeros(Float64, 3))

    if norm(v_r[i]) > 0.0
        signal_angle = dot(vel[i], v_r[i]) / (norm(v_r[i])*norm(vel[i]))

        signal_angle = ifelse( signal_angle < -1, -1, signal_angle)
        signal_angle = ifelse( signal_angle > 1, 1, signal_angle)

        q_r = qrotation(cross(vel[i], v_r[i]), acos(signal_angle) + η * (2.0 * rand() * pi - pi)) * Quaternion(vel[i])

    else

        noise = randn(3)

        signal_angle = dot(noise, p_vel) / (norm(noise)*norm(p_vel))

        signal_angle = ifelse( signal_angle < -1, -1, signal_angle)
        signal_angle = ifelse( signal_angle > 1, 1, signal_angle)

        q_r = qrotation( cross(p_vel, noise), η * acos(signal_angle) ) * Quaternion(p_vel)

    end

	u_vel = normalize([q_r.v1, q_r.v2, q_r.v3])

	vel[i] = u_vel
	pos[i] += u_vel

end
