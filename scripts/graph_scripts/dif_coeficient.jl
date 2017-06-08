### ================================== ###### ================================== ###
### ================================== ###### ================================== ###

using Plots, CollectiveDynamics.DataAnalysis, Polynomials, LaTeXStrings, Quaternions

gr()
pyplot()

gui()
### ================================== ###### ================================== ###

N = 1024
N = 512
N = 256

τ = 7
τ = 6

times = get_times(τ)
x_vals = broadcast(x -> log10(x), times)

### ================================== ###### ================================== ###

folder = "NLOC_DATA"
folder = "NLOC_DATA_3D"
folder = "NLOC_TOP_3D"
folder = "NLOC_MET_3D"
folder = "NLOC_TOP_3D_MEAN"
folder = "TFLOCK_NLOC_DATA"
folder = "TFLOCK_DATA"

folder_path = "$(homedir())/art_DATA/$(folder)/EXP/exp_data_N_$(N)"

### ================================== ###### ================================== ###

make_dir_from_path("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)")
make_dir_from_path("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)/N_$(N)")

eta_folders = readdir(folder_path)

# Para TFLOCK_NLOC_DATA
folder_path = folder_path * "/" * eta_folders[4]

exp_files = filter( x -> ismatch(r"^exp.", x), readdir(folder_path))
nn_files = filter( x -> ismatch(r"^nn_mean.", x), readdir(folder_path))

means = zeros(length(times), length(exp_files))
nn_means = zeros(length(times), length(exp_files))

# std_means = zeros(length(times), length(order_files))

# para NLOC
vals = [ parse(Float64, match(r"^\w+_(\d+\.\d+)", x).captures[1]) for x in exp_files ]

# para TFLOCK
vals = [ parse(Float64, match(r"(\d+\.\d+)\.dat$", x).captures[1]) for x in exp_files ]

### ================================== ###### ================================== ###

for i in sortperm(vals)

    println(i)

    raw_data = reinterpret(Float64, read(folder_path * "/" * exp_files[i]))
    exp_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    raw_data = reinterpret(Float64, read(folder_path * "/" * nn_files[i]))
    nn_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    means[:, i] = mean(exp_data, 2)
    nn_means[:, i] = mean(nn_data, 2)
    # std_means[:, i] = std(exp_data, 2)

end

### ================================== ###### ================================== ###

k_lim = 1.0
lim = 350
m_s = 6

### ================================== ###### ================================== ###

exp_p   = plot(times, hcat([means[:,i] for i in sortperm(vals)]...), lab = [vals[i] for i in sortperm(vals)]', xscale = :log10, yscale = :log10, xlabel = L"t", ylabel = L"\langle r^*(t) \rangle", size = (800,600))

savefig("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)/N_$(N)/r_t.png")

nn_p   = plot(times, hcat([nn_means[:,i] for i in sortperm(vals)]...), lab = [vals[i] for i in sortperm(vals)]', xscale = :log10, yscale = :log10, xlabel = L"t", ylabel = L"\langle r_{nn}^*(t) \rangle", size = (800,600))

savefig("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)/N_$(N)/r_nn_t.png")

### ================================== ###### ================================== ###

k_vals = [vals[i] for i in sortperm(vals)]

r_vals = broadcast(x -> log10(x), means)
r_nn_vals = broadcast(x -> log10(x), nn_means)

r_fit_p   = scatter(x_vals[lim:end], hcat([r_vals[lim:end,i] for i in sortperm(vals)]...), lab = [vals[i] for i in sortperm(vals)]', xlabel = L"t", ylabel = L"\langle r^*(t) \rangle_{\kappa}", alpha = 0.5)

r_nn_fit_p   = scatter(x_vals[lim:end], hcat([r_nn_vals[lim:end,i] for i in sortperm(vals)]...), lab = [vals[i] for i in sortperm(vals)]', xlabel = L"t", ylabel = L"\langle r^*(t) \rangle_{\kappa}", alpha = 0.5)

r_fit_vals = [polyfit(x_vals[lim:end], r_vals[lim:end, i], 1) for i in sortperm(vals)]
r_nn_fit_vals = [polyfit(x_vals[lim:end], r_nn_vals[lim:end, i], 1) for i in sortperm(vals)]

### ================================== ###### ================================== ###

dif_c = plot(k_vals[2:end], [exp10(r_fit_vals[i][0]) for i in 2:length(r_fit_vals)], marker = :o, xscale = :log10, xlabel = L"\kappa", ylabel = L"D", title = "Diffusion Coeffcient", titlefont = font(10), label = L"\langle r \rangle", ms = m_s, size = (800,600), yscale = :log10)

plot!(dif_c, k_vals[2:end], [exp10(r_nn_fit_vals[i][0]) for i in 2:length(r_fit_vals)], marker = :rect, label = L"\langle r_{nn} \rangle", ms = m_s)

savefig("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)/N_$(N)/r_dif_coef.png")

dif_e = plot(k_vals[2:end], [r_fit_vals[i][1] for i in 2:length(r_fit_vals)], marker = :o, xscale = :log10, xlabel = L"\kappa", ylabel = L"\alpha", title = "Diffusion Exponent", titlefont = font(10), label = L"\langle r \rangle", ms = m_s, size = (800,600))

plot!(dif_e, k_vals[2:end], [r_nn_fit_vals[i][1] for i in 2:length(r_fit_vals)], marker = :rect, label = L"\langle r_{nn} \rangle", ms = m_s)

savefig("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)/N_$(N)/r_dif_exp.png")

### ================================== ###### ================================== ###
