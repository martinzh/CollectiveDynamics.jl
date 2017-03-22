### ================================== ###
### ================================== ###

using Plots, CollectiveDynamics.DataAnalysis

gr()
pyplot()

N = 1024
N = 256

τ = 7
τ = 6

times = get_times(τ)

folder_path = "$(homedir())/art_DATA/NLOC_DATA/EXP/exp_data_N_$(N)"

eta_folders = readdir(folder_path)

order_files = filter( x -> ismatch(r"order_k_\d+\.\d+.dat", x), readdir(folder_path))
exp_files = filter( x -> ismatch(r"exp_k_\d+\.\d+.dat", x), readdir(folder_path))

means = zeros(length(times), length(order_files))
orders = zeros(length(times), length(order_files))
std_means = zeros(length(times), length(order_files))

i = 1
for i in 1:length(order_files)

    raw_data = reinterpret(Float64, read(folder_path * "/" * exp_files[i]))
    exp_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    raw_data = reinterpret(Float64, read(folder_path * "/" * order_files[i]))
    order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    means[:, i] = mean(exp_data, 2)
    std_means[:, i] = std(exp_data, 2)
    orders[:, i] = mean(order_data, 2)

end

order_p = plot(times, orders, leg = false, xscale = :log10)
exp_p   = plot(times, means, leg   = false, xscale = :log10, yscale = :log10)
exp_p   = plot(times, means, yerror = std_means, leg   = false, xscale = :log10, yscale = :log10, size = (1024,720))

plot(exp_p, order_p, link = :x, layout = @layout [a ;b])

gui()
