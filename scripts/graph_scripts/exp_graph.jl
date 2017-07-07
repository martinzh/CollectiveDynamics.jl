### ================================== ###
### ================================== ###

using Plots, CollectiveDynamics.DataAnalysis, Polynomials, LaTeXStrings, Quaternions

gr()
pyplot()

gui()
### ================================== ###

N = 1024
N = 512
N = 256
N = 128

N = 8
N = 80
N = 800
N = 2400
N = 4000

τ = 7
τ = 6
τ = 4
τ = 3

times = get_times(τ)
x_vals = broadcast(x -> log10(x), times)

### ================================== ###

function calc_derivs(x_vals, vals)

    all_derivs = zeros(length(times), size(vals, 2))

    # j = 5
    for j in 1:size(vals, 2)
        y_vals = broadcast(x -> log10(x), vals[:, j])

        # ajusta valores
        p_fit = polyfit(x_vals, y_vals, 30)

        # calcula derivada de ajuste
        all_derivs[:, j] = polyval(polyder(p_fit), x_vals)
    end
        return all_derivs
end
### ================================== ###

folder = "NLOC_DATA"
folder = "NLOC_DATA_3D"
folder = "NLOC_TOP_3D"
folder = "NLOC_MET_3D"
folder = "NLOC_VSK_3D"
folder = "NLOC_TOP_3D_MEAN"
folder = "TFLOCK_NLOC_DATA"
folder = "TFLOCK_DATA"
folder = "SVM_GRID_3D"

folder_path = "$(homedir())/art_DATA/$(folder)/EXP/exp_data_N_$(N)"
# folder_path = "$(homedir())/exp_DATA/$(folder)/EXP/exp_data_N_$(N)"

eta_folders = readdir(folder_path)

# Para TFLOCK_NLOC_DATA
folder_path = folder_path * "/" * eta_folders[4]

make_dir_from_path("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)")
make_dir_from_path("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)/N_$(N)")

### ================================== ###### ================================== ###

order_files = filter( x -> ismatch(r"^order.", x), readdir(folder_path))
exp_files = filter( x -> ismatch(r"^exp.", x), readdir(folder_path))
nn_files = filter( x -> ismatch(r"^nn_mean.", x), readdir(folder_path))

means = zeros(length(times), length(order_files))
orders = zeros(length(times), length(order_files))
nn_means = zeros(length(times), length(order_files))

std_means = zeros(length(times), length(order_files))

### ================================== ###### ================================== ###


# para VSK
capt = [match(r"_(\d+\.\d+)\.\w+$|_(\d+\.\d+e-5)\.\w+$", f).captures for f in order_files]
vals = [parse(Float64, vcat(capt...)[i]) for i in find(x -> x != nothing, vcat(capt...))]

# para NLOC
vals = [ parse(Float64, match(r"^\w+_(\d+\.\d+)", x).captures[1]) for x in order_files ]

# para TFLOCK
vals = [ parse(Float64, match(r"(\d+\.\d+)\.dat$", x).captures[1]) for x in order_files ]

### ================================== ###### ================================== ###

for i in sortperm(vals)

    println(i)

    raw_data = reinterpret(Float64, read(folder_path * "/" * exp_files[i]))
    exp_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    raw_data = reinterpret(Float64, read(folder_path * "/" * order_files[i]))
    order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    raw_data = reinterpret(Float64, read(folder_path * "/" * nn_files[i]))
    nn_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    means[:, i] = mean(exp_data, 2)
    std_means[:, i] = std(exp_data, 2)
    orders[:, i] = mean(order_data, 2)
    nn_means[:, i] = mean(nn_data, 2)

end

### ================================== ###### ================================== ###

k_lim = 1.0

y_l = 0.97 #NLOC_DATA
y_h = 1.0 #NLOC_DATA

### ================================== ###### ================================== ###

gui()

# order_p = plot(times, hcat([orders[:,i] for i in sortperm(vals)]...), lab = [vals[i] for i in sortperm(vals)]', xscale = :log10, xlabel = L"t", ylabel = L"\Psi_{\kappa}(t)", size = (600,600), grid = false)

order_p = plot(times, hcat([orders[:,i] for i in sortperm(vals)]...), lab = [vals[i] for i in sortperm(vals)]', xscale = :log10 , grid = false, tickfont = font(12), size = (600,600), aspect_ratio = :equal)

order_p = plot(times, hcat([orders[:,i] for i in sortperm(vals)]...), lab = [vals[i] for i in sortperm(vals)]', xscale = :log10)

savefig("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)/N_$(N)/order_t.png")
### ================================== ###### ================================== ###

psi_plot = plot(vals[sortperm(vals)], [orders[end, sortperm(vals)]], marker = :o,  xlims = (exp10(-5), k_lim), leg = false, xlabel = L"\kappa", ylabel = L"\Psi(\kappa)", xscale = :log10, size = (800,600))
plot!(psi_plot, vals[sortperm(vals)], orders[end, :], marker = :o,  xlims = (k_lim, vals[sortperm(vals)][end]), ylim = (y_l, y_h),  leg = false, inset_subplots = [(1, bbox(0.5w,0.55h,0.45w,0.35h))], subplot=2)

### ================================== ###

psi_plot = plot(vals[sortperm(vals)], [orders[end, sortperm(vals)]], marker = :o,  xlims = (exp10(-2), k_lim), leg = false, grid = false, size = (600,600), tickfont = font(16), aspect_ratio = :equal)

plot!(psi_plot, vals[sortperm(vals)], orders[end, :], marker = :o,  xlims = (k_lim, vals[sortperm(vals)][end]), ylim = (y_l, y_h),  leg = false, inset_subplots = [(1, bbox(0.65w,0.55h,0.25w,0.25h))], subplot=2, grid = false, tickfont = font(10))

savefig("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)/N_$(N)/order_k.png")

N_psi_plot = plot(vals[sortperm(vals)], [orders[end, sortperm(vals)]], marker = :o,  xlims = (exp10(-3), k_lim), xlabel = L"\kappa", ylabel = L"\Psi(\kappa)", xscale = :log10, size = (800,600), lab = "N = $(N)")

plot!(N_psi_plot, vals[sortperm(vals)], [orders[end, sortperm(vals)]], marker = :o,  xlims = (exp10(-3), k_lim), xlabel = L"\kappa", ylabel = L"\Psi(\kappa)", xscale = :log10, size = (800,600), lab = "N = $(N)")

savefig("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)/N_$(N)/order_k_size.png")

### ================================== ###### ================================== ###

exp_p   = plot(times, hcat([means[:,i] for i in sortperm(vals)]...), lab = [vals[i] for i in sortperm(vals)]', xscale = :log10, yscale = :log10, xlabel = L"t", ylabel = L"\langle r^*(t) \rangle", size = (800,600))

exp_p   = plot(times, hcat([means[:,i] for i in sortperm(vals)]...), lab = [vals[i] for i in sortperm(vals)]', xscale = :log10, yscale = :log10, xlims = (exp10(2), exp10(6)), xlabel = L"t", ylabel = L"\langle r^*(t) \rangle_{\kappa}")

# exp_p   = plot(times, means, yerror = std_means, leg   = false, xscale = :log10, yscale = :log10, size = (1024,720))

savefig("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)/N_$(N)/r_t.png")
### ================================== ###### ================================== ###

nn_p   = plot(times, hcat([nn_means[:,i] for i in sortperm(vals)]...), lab = [vals[i] for i in sortperm(vals)]', xscale = :log10, yscale = :log10, xlabel = L"t", ylabel = L"\langle r_{nn}^*(t) \rangle", size = (800,600))

plot!(nn_p, times, fill(10.0, length(times)), linestyle = :dot, linecolor = :red, lab = "r_0")

nn_p   = plot(times, hcat([nn_means[:,i] for i in sortperm(vals)]...), lab = [vals[i] for i in sortperm(vals)]', xscale = :log10, xlims = (exp10(2), exp10(6)), xlabel = L"t", ylabel = L"\langle r_{n.n}(t) \rangle", lw = 1.5, alpha = 0.7)

savefig("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)/N_$(N)/r_nn.png")
### ================================== ###### ================================== ###

exp_d_plot = plot(times, hcat([calc_derivs(x_vals, means)[:, i] for i in sortperm(vals)]...), xscale = :log10, xlims = (exp10(1), exp10(7)), ylims = (0, 1.7), lab = [vals[i] for i in sortperm(vals)]', xlabel = L"t", ylabel = L" \frac{\mathrm{d}\langle r_{ij}(t) \rangle}{\mathrm{d}t} ")

nn_d_plot = plot(times, hcat([calc_derivs(x_vals, nn_means)[:, i] for i in sortperm(vals)]...), xscale = :log10, xlims = (exp10(1), exp10(7)), ylims = (0, 1.2), lab = [vals[i] for i in sortperm(vals)]', xlabel = L"t", ylabel = L" \frac{\mathrm{d}\langle r_{nn}(t) \rangle}{\mathrm{d}t} ", legend = :topleft)

### ================================== ###### ================================== ###

plot(exp_p, order_p, link = :x, layout = @layout [a ;b])

plot(exp_p, nn_p, order_p, layout = @layout [a b c])
plot(exp_p, nn_p, layout = @layout [a b])
plot(exp_p, exp_d_plot, layout = @layout [a b])
plot(nn_p, nn_d_plot, layout = @layout [a b])
plot(order_p, psi_plot, layout = @layout [a b])

gui()
