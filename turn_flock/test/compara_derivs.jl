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

### ================================== ###

function get_times(T)

    times = [convert(Int, exp10(i)) for i in 0:T]

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

    return tau
end

### ================================== ###
function deriv_simple(Y, X)

    N = length(Y)

    δY = zeros(N-1)

    for i in 1:(N-1)

        h = X[i+1] - X[i]
        δY[i] = (Y[i+1] - Y[i]) / h

    end

    return δY
end
### ================================== ###

function deriv_sym(Y, X)

    N = length(Y)

    δY = zeros(N)

    for i in 2:(N-1)

        # h = X[i+1] - X[i-1]
        # δY[i] = (Y[i+1] - Y[i-1]) / (2 * h)

        if i == 1
            h = X[i+1] - X[i]
            δY[i] = (Y[i+1] - Y[i]) / (h)
        elseif i == N
            h = X[i] - X[i-1]
            δY[i] = (Y[i] - Y[i-1]) / (h)
        else
            h = X[i+1] - X[i-1]
            δY[i] = (Y[i+1] - Y[i-1]) / (2 * h)
        end
    end

    return δY
end

### ================================== ###

function deriv_cen(Y, X)

    N = length(Y)

    δY = zeros(N-1)

    for i in 3:(N-2)
        h = X[i+1] - X[i-1]
        δY[i] = (Y[i-2] - 8*Y[i-1] + 8*Y[i+1] - Y[i+2]) / (12 * h)
    end

    return δY
end

### ================================== ###
function deriv_interpol(fX, X)

    n = length(fX)

    δfX = zeros(n)

    for i in 1:3:(n-3)

        h_1 = X[i+1] - X[i]
        h_2 = X[i+2] - X[i+1]

        δfX[i] = ((2*h_1 + h_2)/(h_1*(h_1+h_2))) * fX[i] + ((h_1 + h_2)/(h_1 * h_2)) * fX[i+1] - ((h_1)/(h_2 * (h_1 + h_2))) * fX[i+2]
        δfX[i + 1] = -(h_2/(h_1*(h_1+h_2))) * fX[i] - ((h_1 - h_2)/(h_1 * h_2)) * fX[i+1] - (h_1/(h_2 * (h_1 + h_2))) * fX[i+2]
        δfX[i + 2] = (h_2/(h_1*(h_1+h_2))) * fX[i] - ((h_1 + h_2)/(h_1 * h_2)) * fX[i+1] + ((h_1 + 2*h_2)/(h_2 * (h_1 + h_2))) * fX[i+2]

    end
    return δfX
end

### ================================== ###

using Plots

N = 100
N = 1024
N = 512

T = 4

η = "0.3"

folder_path = "$(homedir())/GitRepos/art_nonLocDiff/cvgn_model/DATA_1/data_N_$(N)/data_N_$(N)_eta_$(η)"
folder_path = "$(homedir())/GitRepos/CollectiveDynamics/cvgn_model/DATA_1/data_N_$(N)/data_N_$(N)_eta_$(η)"

tau = get_times(T)
Δtau = diff(tau)
Δs = unique(Δtau)

### ================================== ###

raw_data = reinterpret(Float64,read(folder_path * "/" * "pos_1.dat"))

pos_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))

pos_data[:, 1:3:3N] = pos_data[:, 1:3:3N] .- mean(pos_data[:, 1:3:3N], 2)
pos_data[:, 2:3:3N] = pos_data[:, 2:3:3N] .- mean(pos_data[:, 2:3:3N], 2)
pos_data[:, 3:3:3N] = pos_data[:, 3:3:3N] .- mean(pos_data[:, 3:3:3N], 2)

tr_pos_data = transpose(pos_data)

means = [mean(calc_rij_vect(tr_pos_data[:, i])) for i in 1:size(tr_pos_data,2)]

deriv = Array{Float64}[]

for j in 1:length(Δs)
    filt_x = [tau[i] for i in findin(Δtau, Δs[j])]
    filt_y = [means[i] for i in findin(Δtau, Δs[j])]

    deriv = vcat(deriv, deriv_simple(filt_y, Δs[j]) )
    # deriv = vcat(deriv, deriv_simple(broadcast(x -> log10(x), filt_y), log10(Δs[j])) )
end

exp = plot(tau, means, markershape = :o, xscale = :log10)
der = plot(tau, deriv, markershape = :o, xscale = :log10)

exp = plot(tau, means, markershape = :o)
der = plot(tau, deriv, markershape = :o)

plot(exp, der, layout = @layout([expan ; der_sym]) )

### ================================== ###

gui()

expan = plot(tau, sorted_all_means,
    xlabel = "tiempo",
    ylabel = "<rij(t)>",
    xscale = :log10,
    yscale = :log10,
    # markershape = :auto,
    # markersize = 6,
    # markeralpha = 0.7,
    label = transpose(sorted_keys),
    legend = :topleft)

der = plot(tau, sorted_all_derivs,
    xlabel = "tiempo",
    ylabel = "<rij(t)>'",
    xscale = :log10,
    markershape = :auto,
    markersize = 6,
    markeralpha = 0.7,
    label = transpose(sorted_keys),
    legend = :topleft)

der_sym = plot(tau, sorted_all_derivs_sym,
    xlabel = "tiempo",
    ylabel = "<rij(t)>'",
    xscale = :log10,
    # markershape = :auto,
    # markersize = 6,
    # markeralpha = 0.7,
    label = transpose(sorted_keys),
    legend = :topleft)

der_cen = plot!(tau, sorted_all_derivs_cen,
    xlabel = "tiempo",
    ylabel = "<rij(t)>'",
    xscale = :log10,
    markershape = :auto,
    markersize = 6,
    markeralpha = 0.7,
    label = transpose(sorted_keys),
    legend = :topleft)

### ================================== ###
plot(expan, deriv, deriv_sym, deriv_cen,
    layout = @layout([expan deriv ; order fase_orden ]),
    size = [1024, 1024])

plot(expan, der_sym, layout = @layout([expan der_sym]) )
# plot(order, fase_orden, layout = @layout([order ; fase_orden]) )
# plot(tau, means[10], xscale = :log10, yscale = :log10)

# savefig("$(homedir())/Google\ Drive/proyecto_martin/imagenes/expansion_N_$(N).png")

### ================================== ###
