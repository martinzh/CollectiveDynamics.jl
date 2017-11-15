
using Distributions

### ============== ### ============== ### ============== ###

addprocs(4)

### ============== ### ============== ### ============== ###

κ   = parse(Float64, ARGS[2]) # average non-local interactions
ω   = parse(Float64, ARGS[3]) # interactions relative weight
T   = parse(Int, ARGS[4]) # integration time steps
rep = parse(Int, ARGS[5])

### ============== ### ============== ### ============== ###

N = 4096
N = 256
ρ = 0.3
η = 0.15
v0 = 1.0
dt = 1.0
l = 0.5

κ = 1.0

L  = cbrt(N / ρ) # size of box

r0 = (v0 * dt) / l # local interaction range

κ_dist = Poisson(κ)

pos = SharedArray{Float64}(3*N) # particles positions
vel = SharedArray{Float64}(3*N) # array of particles' velocities

# k_r = SharedArray{Float64}(N) # local metric interactions
# k_n = SharedArray{Float64}(N) # non local topological interactions

v_r = SharedArray{Float64}(3*N) # local metric interactions
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

# for i in 0:N-1
#     print(i+1,"\t", i*N+1,"\t",(i+1)*N,"\t")
#
#     sh_n = find(x -> x > zero(Float64), Symmetric(R_ij, :L)[i*N+1:(i+1)*N])
#     ln_n = find(x -> x < zero(Float64), Symmetric(R_ij, :L)[i*N+1:(i+1)*N])
#
#     print(length(sh_n))
#     print("\t")
#     println(length(ln_n))
# end

@time @parallel for i in 0:N-1

    # print(i+1,"|\t")

    v_r[3i+1] = 0.0
    v_r[3i+2] = 0.0
    v_r[3i+3] = 0.0

    v_n[3i+1] = 0.0
    v_n[3i+2] = 0.0
    v_n[3i+3] = 0.0

    sh_n = find(x -> x > zero(Float64), Symmetric(R_ij, :L)[i*N+1:(i+1)*N])
    k_sh = length(sh_n)

    ln_n = find(x -> x < zero(Float64), Symmetric(R_ij, :L)[i*N+1:(i+1)*N])

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

@time @parallel for i in 1:length(R_ij)
    R_ij[i] = 0.0
end

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

### ============== ### ============== ### ============== ###
