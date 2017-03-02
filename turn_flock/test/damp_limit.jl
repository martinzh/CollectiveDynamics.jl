### ============== ### ============== ### ============== ###
## Coputation of overdamping limit
## inertial spin model Reference: Cavagna...
## Martin Zumaya Hernandez
## 2 / 3 / 2017
### ============== ### ============== ### ============== ###

using Plots
pyplot()
gr()

### ============== ### ============== ### ============== ###
### COMPUTE RELATIVE DISTANCES MATRIX
### ============== ### ============== ### ============== ###

function calc_nij(pos, nij)

    N = size(nij, 1)

    # compute nij entries
    for i in 1:N, j in (i+1):N
        nij[i,j] = norm(pos[i] - pos[j])
        nij[j,i] = nij[i,j]
    end
end

### ============== ### ============== ### ============== ###

N = 256

ρ   = 0.3
J   = 0.8
χ   = 1.25
v0  = 0.1
n_c = 6

L = cbrt(N / ρ) # 3D

nij = zeros(Float64, N, N)

a = zeros(1000)

for k in 1:1000
    pos = [ [2*rand()*L - L, 2*rand()*L - L, 2*rand()*L - L] for i in 1:N ]
    calc_nij(pos, nij)
    a[k] = mean([mean( [ norm(pos[i] - pos[j]) for j in findin(nij[:,i], sort(nij[:,i])[2:n_c+1]) ] ) for i in 1:N])
end

limit = sqrt(J*n_c*χ)* mean(a) / L

plot(rand(100), leg = :topleft)
