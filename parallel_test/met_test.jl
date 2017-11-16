

### ============== ### ============== ### ============== ###

addprocs(4)

### ============== ### ============== ### ============== ###

@everywhere using Distributions, Quaternions

### ============== ### ============== ### ============== ###

κ   = parse(Float64, ARGS[2]) # average non-local interactions
ω   = parse(Float64, ARGS[3]) # interactions relative weight
T   = parse(Int, ARGS[4]) # integration time steps
rep = parse(Int, ARGS[5])
i_f = ARGS[6]

### ============== ### ============== ### ============== ###

@everywhere N = 4096
@everywhere N = 256

@everywhere κ = 0.1
@everywhere ω = 0.5

@everywhere ρ = 0.3
@everywhere η = 0.15
@everywhere v0 = 1.0
@everywhere dt = 1.0
@everywhere l = 0.5

rep = 1

@everywhere L  = cbrt(N / ρ) # size of box
@everywhere r0 = (v0 * dt) / l # local interaction range

@everywhere κ_dist = Poisson(κ)

@everywhere n_nl = 6

pos = SharedArray{Float64}(3*N) # particles positions
vel = SharedArray{Float64}(3*N) # array of particles' velocities

v_r = SharedArray{Float64}(3*N) # local metric interactions
v_n = SharedArray{Float64}(3*N) # non local topological interactions

R_ij = SharedArray{Float64}(N,N)

# k_r = SharedArray{Float64}(N) # local metric interactions
# k_n = SharedArray{Float64}(N) # non local topological interactions

### ============== ### ============== ### ============== ###
### INITIALIZATION FROM FILE
### ============== ### ============== ### ============== ###

# pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))
# vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

data_path = "$(homedir())/art_DATA/NLOC_P_MET_3D/DATA/data_N_$(N)/data_N_$(N)_k_$(κ)_w_$(ω)"

raw_data = reinterpret(Float64, read(data_path * "/pos_$(rep).dat"))
pos_data = raw_data[(end-3N+1):end]

raw_data = reinterpret(Float64, read(data_path * "/vel_$(rep).dat"))
vel_data = raw_data[(end-3N+1):end]

@time for i in 1:length(pos)
    pos[i] = pos_data[i]
    vel[i] = vel_data[i]
end

pos
vel

### ============== ### ============== ### ============== ###
### INITIALIZATION
### ============== ### ============== ### ============== ###

@time for i in 1:length(pos)
    pos[i] = 2*rand()*L - L
    vel[i] = 2*rand() - 1
end

@time for i in 1:3:length(vel)
    norm = sqrt(vel[i]^2 + vel[i+1]^2 + vel[i+2]^2)
    vel[i] /= norm
    vel[i+1] /= norm
    vel[i+2] /= norm
end

### ============== ### ============== ### ============== ###
### COMPUTE RELATIVE DISTANCES WITH METRIC THRESHOLD
### ============== ### ============== ### ============== ###

@time @parallel for i in 1:3:length(pos)
        for j in (i+3):3:length(pos)
            d = (pos[i]-pos[j])^2 + (pos[i+1]-pos[j+1])^2 + (pos[i+2]-pos[j+2])^2
            d <= r0^2 && d > 0.0 ? R_ij[div(j,3)+1,div(i,3)+1] = 1.0 : R_ij[div(j,3)+1,div(i,3)+1] = -1.0
        end
end

### ============== ### ============== ### ============== ###
### COMPUTE RAW RELATIVE DISTANCES
### ============== ### ============== ### ============== ###

@time @parallel for i in 1:3:length(pos)
        for j in (i+3):3:length(pos)
            R_ij[div(j,3)+1,div(i,3)+1] = (pos[i]-pos[j])^2 + (pos[i+1]-pos[j+1])^2 + (pos[i+2]-pos[j+2])^2
        end
end

R_ij
Symmetric(R_ij, :L)

### ============== ### ============== ### ============== ###
### COMPUTE SHORT AND LONG RANGE METRIC INTERACTIONS
### ============== ### ============== ### ============== ###

@everywhere function compute_interactions(vel::SharedArray,v_r::SharedArray,v_n::SharedArray)

    for id in first(localindexes(vel)):3:last(localindexes(vel))

        i = div(id, 3)
        # print(i+1,"|\t")

        v_r[3i+1] = 0.0
        v_r[3i+2] = 0.0
        v_r[3i+3] = 0.0

        v_n[3i+1] = 0.0
        v_n[3i+2] = 0.0
        v_n[3i+3] = 0.0

        sh_n = find(x -> x > zero(Float64), Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])
        k_sh = length(sh_n)
        # k_sh = sum(sh_n)

        ln_n = find(x -> x < zero(Float64), Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])

        # short-range
        if isempty(sh_n) == false
            for j in sh_n
                # print(j,"\t")
                v_r[3i+1] += vel[3(j-1)+1] / k_sh
                v_r[3i+2] += vel[3(j-1)+2] / k_sh
                v_r[3i+3] += vel[3(j-1)+3] / k_sh
            end
        end

        # println()
        k_ln = rand(κ_dist)
        # println(k_ln)

        # possible long range
        # if isempty(ln_n) == false && k_ln != 0.0
        if k_ln != 0.0
            for j in rand(ln_n, 3)
                # print(j,"\t")
                v_n[3i+1] += vel[3(j-1)+1] / k_ln
                v_n[3i+2] += vel[3(j-1)+2] / k_ln
                v_n[3i+3] += vel[3(j-1)+3] / k_ln
            end
        end

        # println()

    end
end

### ============== ### ============== ### ============== ###
### COMPUTE SHORT AND LONG RANGE TOPOLOGICAL INTERACTIONS
### ============== ### ============== ### ============== ###

@everywhere function compute_interactions(vel::SharedArray,v_r::SharedArray,v_n::SharedArray)

    for id in first(localindexes(vel)):3:last(localindexes(vel))

        i = div(id, 3)
        # print(i+1,"|\t")

        v_r[3i+1] = 0.0
        v_r[3i+2] = 0.0
        v_r[3i+3] = 0.0

        v_n[3i+1] = 0.0
        v_n[3i+2] = 0.0
        v_n[3i+3] = 0.0


        # first neighbors
        sh_n = findin(Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N], sort(Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])[2:k_sh+1])

        # next neighbors
        ln_n = findin(Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N], sort(Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])[k_sh+2:end])

        # short-range
        for j in sh_n
            # print(j,"\t")
            v_r[3i+1] += vel[3(j-1)+1] / k_sh
            v_r[3i+2] += vel[3(j-1)+2] / k_sh
            v_r[3i+3] += vel[3(j-1)+3] / k_sh
        end

        # println()
        k_ln = rand(κ_dist)
        # println(k_ln)

        # possible long range
        # if isempty(ln_n) == false && k_ln != 0.0
        if k_ln != 0.0
            for j in rand(ln_n, k_ln)
                # print(j,"\t")
                v_n[3i+1] += vel[3(j-1)+1] / k_ln
                v_n[3i+2] += vel[3(j-1)+2] / k_ln
                v_n[3i+3] += vel[3(j-1)+3] / k_ln
            end
        end

        # println()

    end
end

### ============== ### ============== ### ============== ###
### UPDATE PARTICLE'S POSITIONS AND VELOCITIES
### ============== ### ============== ### ============== ###

@everywhere function update_particles(pos::SharedArray, vel::SharedArray,v_r::SharedArray,v_n::SharedArray)

    for id in first(localindexes(vel)):3:last(localindexes(vel))

        i = div(id, 3)
        # print(i+1,"|\t")

        signal = ω * [v_r[3i+1] , v_r[3i+2], v_r[3i+3]] + (1.0 - ω) * [v_n[3i+1] , v_n[3i+2], v_n[3i+3]]
        signal_angle = 0.0

        p_vel = [vel[3i+1] , vel[3i+2], vel[3i+3]]

        norm(signal) != zero(Float64) ? signal_angle = acos(dot(p_vel, signal) / norm(signal)) : signal_angle = 0.0

        q_r = Quaternion(zeros(Float64, 3))

        if norm(signal) != zero(Float64)
            q_r = qrotation(cross(p_vel, signal), signal_angle + η * (2.0 * rand() * pi - pi)) * Quaternion(p_vel)
        else
            noise = randn(3)
            q_r = qrotation(cross(p_vel, noise), η * acos(dot(normalize(noise), p_vel)) ) * Quaternion(p_vel)
        end

        u_vel = normalize([q_r.v1, q_r.v2, q_r.v3])

        vel[3i+1] = u_vel[1]
        vel[3i+2] = u_vel[2]
        vel[3i+3] = u_vel[3]

        pos[3i+1] += u_vel[1]
        pos[3i+2] += u_vel[2]
        pos[3i+3] += u_vel[3]

    end
end

### ============== ### ============== ### ============== ###
### SYSTEM UPDATE
### ============== ### ============== ### ============== ###

function evolve_system()

    @parallel for i in 1:3:length(pos)
            for j in (i+3):3:length(pos)
                R_ij[div(j,3)+1,div(i,3)+1] = 0.0
                d = (pos[i]-pos[j])^2 + (pos[i+1]-pos[j+1])^2 + (pos[i+2]-pos[j+2])^2
                d <= r0^2 && d > 0.0 ? R_ij[div(j,3)+1,div(i,3)+1] = 1.0 : R_ij[div(j,3)+1,div(i,3)+1] = -1.0
            end
    end

    for p in workers()
        remotecall_fetch(compute_interactions, p, vel, v_r, v_n)
    end

    for p in workers()
        remotecall_fetch(update_particles, p, pos, vel, v_r, v_n)
    end
end

@time evolve_system()

### ============== ### ============== ### ============== ###
### TEST ZONE
### ============== ### ============== ### ============== ###

for i in 1:N
    k_r[i] = 0.0
end

for i in 1:length(v_r)
    v_r[i] = 0.0
    v_n[i] = 0.0
end

for i in 1:length(R_ij)
    R_ij[i] = 0.0
end

[norm([v_r[i],v_r[i+1],v_r[i+2]]) for i in 1:3:length(v_r)]
[norm([v_n[i],v_n[i+1],v_n[i+2]]) for i in 1:3:length(v_r)]
[norm([vel[i],vel[i+1],vel[i+2]]) for i in 1:3:length(v_r)]

### ============== ### ============== ### ============== ###

for i in 1:3:192
    println(i,"\t", div(i,3),"\t", (div(i,3)*N)+1, "\t", (div(i,3)+1)*N)
end

for i in procs()
    println(remotecall_fetch(localindexes, i, vel))
    # println(remotecall_fetch(localindexes, i, v_r))
    # println(remotecall_fetch(localindexes, i, v_n))
end

indexpids(pos)

### ============== ### ============== ### ============== ###

@time for p in workers()
    remotecall_fetch(compute_interactions, p, vel, v_r, v_n)
end

@time for p in workers()
    remotecall_fetch(update_particles, p, pos, vel, v_r, v_n)
end

v_r
v_n
pos
vel

### ============== ### ============== ### ============== ###

i = 1

# first neighbors
findin(Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N], sort(Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])[2:n_nl+1])

# next neighbors
rand(findin(Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N], sort(Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])[n_nl+2:end]), 3)

### ============== ### ============== ### ============== ###

Ti = 0
Tf = 5

times = [convert(Int, exp10(i)) for i in Ti:Tf]

tes = []

for i in 1:(length(times) - 1)

    if i > 1

        println(i)

        for t in (times[i]+1):times[i+1]

            if t % times[i] == 0 || t % times[i-1] == 0
                # println("//////// ", t)
                push!(tes, t)
            end
        end

    else

        for t in (times[i]+1):times[i+1]

            if t % times[i] == 0
                # println("//////// ", t)
                push!(tes, t)
            end
        end

    end

end

length(tes)
