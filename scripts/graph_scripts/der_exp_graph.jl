
using Plots, CollectiveDynamics.DataAnalysis, Polynomials, LaTeXStrings, Quaternions
### ================================== ###

gui()
gr()
pyplot()

N = 1024
N = 256

τ = 7
τ = 6

times = get_times(τ)

### ================================== ###

folder = "NLOC_DATA"
folder = "NLOC_DATA_3D"
folder = "NLOC_TOP_3D"
folder = "NLOC_TOP_3D_MEAN"
folder = "TFLOCK_NLOC_DATA"

folder_path = "$(homedir())/art_DATA/$(folder)/EXP/exp_data_N_$(N)"

eta_folders = readdir(folder_path)

order_files = filter( x -> ismatch(r"^order.", x), readdir(folder_path))
exp_files = filter( x -> ismatch(r"^exp.", x), readdir(folder_path))
nn_files = filter( x -> ismatch(r"^nn_mean.", x), readdir(folder_path))
vel_files = filter( x -> ismatch(r"^vel_mean.", x), readdir(folder_path))

# para TFLOCK_NLOC_DATA
order_files = filter( x -> ismatch(r"^order.", x), readdir(folder_path * "/" * eta_folders[1]))
exp_files = filter( x -> ismatch(r"^exp.", x), readdir(folder_path * "/" * eta_folders[1]))
nn_files = filter( x -> ismatch(r"^nn_mean.", x), readdir(folder_path * "/" * eta_folders[1]))
vel_files = filter( x -> ismatch(r"^vel_mean.", x), readdir(folder_path * "/" * eta_folders[1]))

# k_vals = [ match(r"(\d+\.\d+)\w+\d+\.\d+.dat$", x).captures[1] for x in order_files]

# para TFLOCK_NLOC_DATA
k_vals = [ parse(Float64,match(r"(\d+\.\d+)\.dat$", x).captures[1]) for x in order_files]

# para NLOC_DATA
k_vals = [ parse(Float64, match(r"(\d+\.\d+).dat$", x).captures[1]) for x in order_files]

# para  NLOC_TOP
k_vals = [ parse(Float64,match(r"^\w+(\d+\.\d+)_.", x).captures[1]) for x in order_files]

### ================================== ###

means = zeros(length(times), length(order_files))
nn_means = zeros(length(times), length(order_files))
vel_means = zeros(2*length(times), length(order_files))
orders = zeros(length(times), length(order_files))
std_means = zeros(length(times), length(order_files))

i = 1
k = 1
for i in sortperm(k_vals)

    raw_data = reinterpret(Float64, read(folder_path * "/" * exp_files[i]))

    # para TFLOCK_NLOC_DATA
    # raw_data = reinterpret(Float64, read(folder_path * "/" * eta_folders[1] * "/" * exp_files[i]))

    exp_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    raw_data = reinterpret(Float64, read(folder_path * "/" * order_files[i]))

    # para TFLOCK_NLOC_DATA
    # raw_data = reinterpret(Float64, read(folder_path * "/" * eta_folders[1] * "/" * order_files[i]))

    order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    raw_data = reinterpret(Float64, read(folder_path * "/" * nn_files[i]))

    nn_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    # raw_data = reinterpret(Float64, read(folder_path * "/" * vel_files[i]))
    #
    # vel_data = reshape(raw_data, 2*length(times), div(length(raw_data), 2*length(times)))

    means[:, k] = mean(exp_data, 2)
    std_means[:, k] = std(exp_data, 2)
    orders[:, k] = mean(order_data, 2)
    nn_means[:, k] = mean(nn_data, 2)
    # vel[:, k] = mean(order_data, 2)

    k += 1
end
### ================================== ###

all_derivs = zeros(length(times), size(means, 2))

x_vals = broadcast(x -> log10(x), times)

# j = 5
for j in 1:size(means, 2)
    y_vals = broadcast(x -> log10(x), means[:, j])

    # ajusta valores
    p_fit = polyfit(x_vals, y_vals, 30)

    # calcula derivada de ajuste
    all_derivs[:, j] = polyval(polyder(p_fit), x_vals)
end

### ================================== ###

s_lag = 10 #NLOC_DATA
s_lag = 3 #NLOC_TOP_3D
s_lag = 3 #NLOC_TOP_3D

e_lag = 150 #NLOC_DATA
e_lag = 0 #NLOC_DATA
e_lag = 0 #NLOC_DATA

k_lim = 8 #NLOC_DATA
k_lim = 0.35 #NLOC_TOP_3D
k_lim = exp10(0) #TFLOCK_NLOC_DATA

y_l = 0.925 #NLOC_DATA
y_h = 0.975 #NLOC_DATA

y_l = 0.975 #NLOC_TOP_3D
y_h = 0.99 #NLOC_TOP_3D

y_l = 0.99 #TFLOCK_NLOC_DATA
y_h = 1.01 #TFLOCK_NLOC_DATA

order_p = plot(times, orders, lab = reshape(k_vals[sortperm(k_vals)], 1, length(k_vals)), xscale = :log10, leg = false, xlabel = L"t", ylabel = L"\Psi_{\kappa}(t)")

exp_p   = plot(times, means, lab = reshape(k_vals[sortperm(k_vals)], 1, length(k_vals)), xscale = :log10, yscale = :log10, leg = false, xlabel = L"t", ylabel = L"\langle r_{ij}(t) \rangle")

nn_p   = plot(times, nn_means, lab = reshape(k_vals[sortperm(k_vals)], 1, length(k_vals)), xscale = :log10, yscale = :log10, leg = false, xlabel = L"t", ylabel = L"\langle r_{nn}(t) \rangle")


plot(exp_p, nn_p, layout = @layout [a b])

exp_p   = plot(times, means, lab = reshape(k_vals[sortperm(k_vals)], 1, length(k_vals)), xscale = :log10, yscale = :log10, leg = false, xlabel = L"t", ylabel = L"\langle r_{ij} \rangle", xlims = (5exp10(2), exp10(6)))


exp_p   = plot(times, means, yerror = std_means, leg = false, xscale = :log10, yscale = :log10, size = (1024,720))

# d_p = plot(broadcast(x -> exp10(x), x_vals[s_lag:end-e_lag]), all_derivs[s_lag:end-e_lag, :], lab = reshape(k_vals[sortperm(k_vals)], 1, length(k_vals)), ylims = (0, 1.7), xscale = :log10, leg = false)
d_p = plot(broadcast(x -> exp10(x), x_vals[s_lag:end]), all_derivs[s_lag:end, :], lab = reshape(k_vals[sortperm(k_vals)], 1, length(k_vals)), ylims = (0, 1.7), xscale = :log10, leg = false, xlabel = L"t", ylabel = L"\frac{\mathrm{d} \langle r_{ij} \rangle}{\mathrm{d}t}")

d_p = plot(broadcast(x -> exp10(x), x_vals[s_lag:end-e_lag]), all_derivs[s_lag:end-e_lag, :], lab = reshape(k_vals[sortperm(k_vals)], 1, length(k_vals)), ylims = (0, 1.8), xscale = :log10, leg = false, xlabel = L"t", ylabel = L"\frac{\mathrm{d} \langle r_{ij} \rangle}{\mathrm{d}t}", xlims = (5exp10(2), exp10(6)))

taus = [times[findmax(all_derivs[s_lag:end-e_lag,i])[2]] for i in 1:size(all_derivs, 2)]

t_plot = plot(k_vals[sortperm(k_vals)], taus, marker = :o, yscale = :log10,  xlims = (0, k_lim), ylims = (exp10(2), 4exp10(4)), leg = false, xlabel = L"\kappa", ylabel = L"\tau")
plot!(t_plot, k_vals[sortperm(k_vals)], taus, marker = :o, yscale = :log10,  xlims = (k_lim, k_vals[sortperm(k_vals)][end]), ylim = (2.5exp10(4), 3.5exp10(4)), leg = false, inset_subplots = [(1, bbox(0.5w,0.55h,0.45w,0.35h))], subplot=2)

psi_plot = plot(k_vals[sortperm(k_vals)], orders[end, :], marker = :o,  xlims = (exp10(-4), k_lim), leg = false, xlabel = L"\kappa", ylabel = L"\Psi(\kappa)", xscale = :log10)
plot!(psi_plot, k_vals[sortperm(k_vals)], orders[end, :], marker = :o,  xlims = (k_lim, k_vals[sortperm(k_vals)][end]), ylim = (y_l, y_h),  leg = false, inset_subplots = [(1, bbox(0.5w,0.55h,0.45w,0.35h))], subplot=2)

plot(exp_p, d_p, order_p, t_plot, psi_plot, layout = @layout [a b c; d e])
plot(exp_p, d_p, order_p, psi_plot, layout = @layout [a b ; d e])

findmax(all_derivs[s_lag:end-e_lag, 1])
