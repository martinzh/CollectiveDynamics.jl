
using Plots, CollectiveDynamics.DataAnalysis
using Polynomials

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


gr()
pyplot()

N = 1024
N = 256

τ = 7
τ = 6

times = get_times(τ)

folder = "NLOC_DATA"
folder = "NLOC_TOP_3D"
folder = "NLOC_TOP_3D_MEAN"

folder_path = "$(homedir())/art_DATA/$(folder)/EXP/exp_data_N_$(N)"

eta_folders = readdir(folder_path)

order_files = filter( x -> ismatch(r"^order.", x), readdir(folder_path))
exp_files = filter( x -> ismatch(r"^exp.", x), readdir(folder_path))

# k_vals = [ match(r"(\d+\.\d+)\w+\d+\.\d+.dat$", x).captures[1] for x in order_files]
k_vals = [ match(r"(\d+\.\d+).dat$", x).captures[1] for x in order_files]

means = zeros(length(times), length(order_files))
orders = zeros(length(times), length(order_files))
std_means = zeros(length(times), length(order_files))

# i = 1
for i in 1:length(order_files)

    raw_data = reinterpret(Float64, read(folder_path * "/" * exp_files[i]))
    exp_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    raw_data = reinterpret(Float64, read(folder_path * "/" * order_files[i]))
    order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    means[:, i] = mean(exp_data, 2)
    std_means[:, i] = std(exp_data, 2)
    orders[:, i] = mean(order_data, 2)

end

order_p = plot(times, orders, lab = k_vals', xscale = :log10)
exp_p   = plot(times, means, lab = k_vals', xscale = :log10, yscale = :log10)
exp_p   = plot(times, means, yerror = std_means, leg   = false, xscale = :log10, yscale = :log10, size = (1024,720))

plot(exp_p, order_p, link = :x, layout = @layout [a b])

gui()

plot(times, means[:,end], xscale = :log10, yscale = :log10)

Δtau = diff(times)
Δs = unique(Δtau)

all_derivs = zeros(length(Δtau), size(means, 2))
all_derivs = zeros(542, size(means, 2))

for k in 1:size(means, 2)

    deriv = Array{Float64}[]

    for j in 1:length(Δs)
        filt_x = [times[i] for i in findin(Δtau, Δs[j])]
        filt_y = [means[i, k] for i in findin(Δtau, Δs[j])]

        # deriv = vcat(deriv, deriv_simple(filt_y, filt_x) )
        deriv = vcat(deriv, deriv_simple(broadcast(x -> log10(x), filt_y), broadcast(x -> log10(x), filt_x)) )
    end

    all_derivs[:, k] = deriv

    end

derivs_p = plot(times, all_derivs, xscale = :log10, lab = k_vals')

plot(exp_p, derivs_p, link = :x, layout = @layout [a b])

gui()

scatter(times, all_derivs[:,1], xscale = :log10, alpha = 0.3)
