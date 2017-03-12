function calc_rij(pos, rij, r0)

    # compute rij entries
    for i in 1:size(rij,1), j in (i+1):size(rij,1)

        d = sqrt((pos[i][1] - pos[j][1])^2 + (pos[i][2] - pos[j][2])^2)

        d < r0 && d > zero(Float64) ? rij[j,i] = one(Float64) : rij[j,i] = zero(Float64)

        rij[i,j] = rij[j,i]
    end

end

function part_inter_vel(vel, Nij)
    # sum(Nij) > zero(Float64) ? sum(map(*, vel, Nij)) / sum(Nij) : sum(map(*, vel, Nij))
    sum(Nij) > zero(Float64) ? sum(map(*, vel, Nij)) / sum(Nij) : zeros(2)
end

function calc_interactions(vel, v_r, v_n, rij, Nij)

    for i in 1:size(Nij, 1)
        v_r[i] = part_inter_vel(vel, view(rij, :, i))
        v_n[i] = part_inter_vel(vel, view(Nij, :, i))
    end

end

function test_Nij(vel, v_n, Nij)

    for i in 1:size(Nij, 1)
        v_n[i] = part_inter_vel(vel, view(Nij, :, i))
    end

end

function test_sparse_Nij(vel, v_n, sp_Nij)

    rows = rowvals(sp_Nij)

    for i in 1:size(sp_Nij, 1)

        range = nzrange(sp_Nij, i)

        length(range) > 0 ? v_n[i] = mean([vel[rows[j]] for j in range ] ) : v_n[i] = zeros(2)

    end

end


function set_Nij(p, Nij)

    N = size(Nij, 1)

    for i in 1:N, j in union(1:i-1, i+1:N)
        rand() < p ? Nij[j, i] = 1 : Nij[j, i] = 0
    end
end

N = 10
κ = 2

p = κ / (N-1)

ρ  = 0.5 # density
L  = sqrt(N / ρ) # size of box

v0 = 1.0

# array of random initial particles' postitions
pos = [ [2*rand()*L - L, 2*rand()*L - L] for i in 1:N ]

# array of particles' velocities
vel = v0 * [ normalize([2*rand() - 1, 2*rand() - 1]) for i in 1:N ]

# local metric interactions
v_r = [zeros(2) for i in 1:N]

# non local topological interactions
v_n = [zeros(2) for i in 1:N]

rij = zeros(N, N)
Nij = zeros(N, N)

set_Nij(p, Nij)

sp_Nij = sparse(Nij)
sp_rij = sparse(rij)

@time test_Nij(vel, v_n, Nij)
@time test_Nij(vel, v_r, rij)


@time test_sparse_Nij(vel, v_n, sp_Nij) ######
@time test_sparse_Nij(vel, v_r, sp_rij) ######

@time calc_rij(pos, rij, 10.0)
