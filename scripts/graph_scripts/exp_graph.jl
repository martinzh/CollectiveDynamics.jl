### ================================== ###
### ================================== ###

using Plots, CollectiveDynamics.DataAnalysis, LaTeXStrings

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
folder = "TFLOCK_NLOC_DATA"

folder_path = "$(homedir())/art_DATA/$(folder)/EXP/exp_data_N_$(N)"

eta_folders = readdir(folder_path)

folder_path = folder_path * "/" * eta_folders[1]

order_files = filter( x -> ismatch(r"^order.", x), readdir(folder_path))
exp_files = filter( x -> ismatch(r"^exp.", x), readdir(folder_path))
nn_files = filter( x -> ismatch(r"^nn_mean.", x), readdir(folder_path))

means = zeros(length(times), length(order_files))
orders = zeros(length(times), length(order_files))
nn_means = zeros(length(times), length(order_files))

std_means = zeros(length(times), length(order_files))

i = 1
for i in 1:length(order_files)

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

order_p = plot(times, orders, leg = false, xscale = :log10, xlabel = L"t", ylabel = L"\Psi_{\kappa}(t)")
exp_p   = plot(times, means, leg   = false, xscale = :log10, yscale = :log10, xlabel = L"t", ylabel = L"\langle r_{ij}(t) \rangle")
exp_p   = plot(times, means, yerror = std_means, leg   = false, xscale = :log10, yscale = :log10, size = (1024,720))

nn_p   = plot(times, nn_means, leg   = false, xscale = :log10, yscale = :log10, xlabel = L"t", ylabel = L"\langle r_{n.n}(t) \rangle")

plot(exp_p, order_p, link = :x, layout = @layout [a ;b])

plot(exp_p, nn_p, order_p, layout = @layout [a b c])
plot(exp_p, nn_p, layout = @layout [a b])

gui()
