

### ============== ### ============== ### ============== ###

addprocs(4)

### ============== ### ============== ### ============== ###

κ   = parse(Float64, ARGS[2]) # average non-local interactions
ω   = parse(Float64, ARGS[3]) # interactions relative weight
T   = parse(Int, ARGS[4]) # integration time steps
rep = parse(Int, ARGS[5])

### ============== ### ============== ### ============== ###

N = 4096
ρ = 0.3
η = 0.15
v0 = 1.0
dt = 1.0
l = 0.5

L  = cbrt(N / ρ) # size of box

r0 = (v0 * dt) / l # local interaction range

# κ_dist = Poisson(κ)

pos = SharedArray{Float64}(3*N) # particles positions
vel = SharedArray{Float64}(3*N) # array of particles' velocities

k_r = SharedArray{Float64}(N) # local metric interactions
v_r = SharedArray{Float64}(3*N) # local metric interactions

k_n = SharedArray{Float64}(N) # non local topological interactions
v_n = SharedArray{Float64}(3*N) # non local topological interactions

R_ij = SharedArray{Float64}(N,N)

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
### COMPUTE RELATIVE DISTANCES
### ============== ### ============== ### ============== ###

@time @parallel for i in 1:3:length(pos)
        for j in (i+3):3:length(pos)
            R_ij[div(j,3)+1,div(i,3)+1] = 0.0
            d = (pos[i]-pos[j])^2 + (pos[i+1]-pos[j+1])^2 + (pos[i+2]-pos[j+2])^2
            d <= r0^2 && d > 0.0 ? R_ij[div(j,3)+1,div(i,3)+1] = 1.0 : R_ij[div(j,3)+1,div(i,3)+1] = -1.0
        end
end

R_ij

### ============== ### ============== ### ============== ###
### COMPUTE SHORT AND LONG RANGE INTERACTIONS
### ============== ### ============== ### ============== ###
v_r
v_n

i = 3

@time for i in 0:N-3

    v_r[3i+1] = 0.0
    v_r[3i+2] = 0.0
    v_r[3i+3] = 0.0

    v_n[3i+1] = 0.0
    v_n[3i+2] = 0.0
    v_n[3i+3] = 0.0

    # short-range
    if isempty(find(x -> x > zero(Float64), Symmetric(R_ij, :L)[i*N+1:(i+1)*N])) == false
        for j in find(x -> x > zero(Float64), Symmetric(R_ij, :L)[i*N+1:(i+1)*N])
            v_r[3i+1] += vel[3j+1]
            v_r[3i+2] += vel[3j+2]
            v_r[3i+3] += vel[3j+3]
        end
    end

    # possible long range
    # find(x -> x < zero(Float64), Symmetric(R_ij, :L)[i*N + 1:(i+1)*N])
    if isempty(find(x -> x < zero(Float64), Symmetric(R_ij, :L)[i*N+1:(i+1)*N])) == false
        for j in rand(find(x -> x < zero(Float64), Symmetric(R_ij, :L)[i*N+1:(i+1)*N]), 3)
            v_n[3i+1] += vel[3j+1]
            v_n[3i+2] += vel[3j+2]
            v_n[3i+3] += vel[3j+3]
        end
    end

end

@time @parallel for i in 1:length(R_ij)
    R_ij[i] = 0.0
end


### ============== ### ============== ### ============== ###


function short_range(i)

    k_r[i+1] = 0.0
    v_r[3i+1] = 0.0
    v_r[3i+2] = 0.0
    v_r[3i+3] = 0.0

    for j in 1:3:length(vel)
        d = (pos[3i+1]-pos[j])^2 + (pos[3i+2]-pos[j+1])^2 + (pos[3i+3]-pos[j+2])^2
        if d <= r0^2 && d > 0.0
            k_r[i+1] += 1.0
            v_r[3i+1] += vel[j]
            v_r[3i+2] += vel[j+1]
            v_r[3i+3] += vel[j+2]
        end
    end

    end
end

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

k_r
v_r
vel

@time pmap(short_range, collect(0:N-1))
