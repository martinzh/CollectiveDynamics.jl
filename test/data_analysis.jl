ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"

using StatPlots

N = 1024
κ = 15
T = 7

folder_path = "data_N_$(N)_k_$(κ)"

means = Array{Float64}[]

# r = 5

for r in 1:15

    println(r)

    raw_data = reinterpret(Float64,read(folder_path * "/pos_$(r).dat"))

    pos_data = reshape(raw_data, 2N, div(length(raw_data), 2N))

    push!(means, [mean(calc_rij_vect(pos_data[:, i])) for i in 1:size(pos_data,2)])

end

# raw_data = reinterpret(Float64,read(folder_path * "/vel_$(r).dat"))
#
# vel_data = reshape(raw_data, 2N, div(length(raw_data), 2N))
#
# raw_data = reinterpret(Float64,read(folder_path * "/net_$(r).dat"))
#
# int_net = reshape(raw_data, N, N)
#
# anim = Animation()
#
# for t in 1:size(pos_data, 2)
#
#     println(t)
#
#     x = [pos_data[i, t] for i in 1:2:2N]
#     y = [pos_data[i+1, t] for i in 1:2:2N]
#
#     pts = vec(P2[(x[i],y[i]) for i in 1:N])
#
#     # overlapping graphs
#     # quiver!(pts, quiver = ([vel[i][1] for i in 1:N], [vel[i][2] for i in 1:N]))
#
#     # overwriting graphs
#     quiver(pts, quiver = ([vel_data[i, t] for i in 1:2:2N], [vel_data[i+1, t] for i in 1:2:2N]))
#
#     frame(anim)
#
# end
#
#
# gif(anim, "test.gif", fps = 24)

### ================================== ###
### ================================== ###

function calc_rij_vect(vect)

    n = div(length(vect),2)

    # println(size(vect))

    vec_rij = zeros(Float64, div(n*(n-1),2) )

    k = 1

    for i in 1:2:2n, j in (i+2):2:2n

        vec_rij[k] = sqrt( (vect[i] - vect[j])^2 + (vect[i+1] - vect[j+1])^2)

        k += 1
    end

    # println("Finish calc rij")

    return vec_rij

end

times = [convert(Int, exp10(i)) for i in 0:T]

mean_rij = [mean(calc_rij_vect(pos_data[:, i])) for i in 1:size(pos_data,2)]

tau = Int64[]

for i in 1:(length(times) - 1)

    if i > 1

        for t in (times[i]+1):times[i+1]

            if t % times[i] == 0 || t % times[i-1] == 0
                push!(tau, t)
            end
        end

    else

        for t in (times[i]+1):times[i+1]

            if t % times[i] == 0
                push!(tau, t)
            end
        end

    end

end


plot(tau, mean(hcat(means...), 2), xscale = :log, yscale = :log)
