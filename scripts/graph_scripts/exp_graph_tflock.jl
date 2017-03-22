### ================================== ###### ================================== ###

using Plots, CollectiveDynamics.DataAnalysis

# gr()
pyplot()

# N = 1024
N = 256

# τ = 7
τ = 6

times = get_times(τ)

folder_path = "$(homedir())/art_DATA/TFLOCK_DATA/EXP/exp_data_N_$(N)"
folder_path = "$(homedir())/art_DATA/TFLOCK_NLOC_DATA/EXP/exp_data_N_$(N)"

eta_folders = readdir(folder_path)

### ================================== ###### ================================== ###

f = 1
for f in 1:length(eta_folders)

    println(eta_folders[f])

    order_files = filter( x -> ismatch(r"^order.", x), readdir(folder_path * "/" * eta_folders[f]))
    exp_files = filter( x -> ismatch(r"^exp.", x), readdir(folder_path * "/" * eta_folders[f]))

    means = zeros(length(times), length(order_files))
    orders = zeros(length(times), length(order_files))
    std_means = zeros(length(times), length(order_files))

    noise_vals = [match(r"(\d+\.\d+)\.dat$", x).captures[1] for x in order_files]

    # i = 1
    for i in 1:length(order_files)

        raw_data = reinterpret(Float64, read(folder_path * "/" * eta_folders[f] * "/" * exp_files[i]))
        exp_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

        raw_data = reinterpret(Float64, read(folder_path * "/" * eta_folders[f] * "/" * order_files[i]))
        order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

        means[:, i] = mean(exp_data, 2)
        std_means[:, i] = std(exp_data, 2)
        orders[:, i] = mean(order_data, 2)

    end

    # order_p = plot(times, orders, leg = false, xscale = :log10)
    # exp_p   = plot(times, means, leg   = false, xscale = :log10, yscale = :log10)
    # exp_p   = plot(times, means, yerror = std_means, leg   = false, xscale = :log10, yscale = :log10, size = (1024,720))

    order_p = plot(times, orders, lab = noise_vals', xscale = :log10)
    exp_p   = plot(times, means, lab = noise_vals', xscale = :log10, yscale = :log10)
    # exp_p   = plot(times, means, yerror = std_means, xscale = :log10, yscale = :log10, size = (1024,720))

    plot(exp_p, order_p, layout = @layout([a b]), size = [1024, 720])
    # plot(exp_p, order_p, layout = @layout[exp_p order_p])

    savefig("$(homedir())/Google\ Drive/proyecto_martin/imagenes/modelo_cvgn/expansion_N_$(N)_$(eta_folders[f]).png")


end

### ================================== ###### ================================== ###

for f in 1:length(eta_folders)

    println(eta_folders[f])

    order_files = filter( x -> ismatch(r"^order.", x), readdir(folder_path * "/" * eta_folders[f]))
    exp_files = filter( x -> ismatch(r"^exp.", x), readdir(folder_path * "/" * eta_folders[f]))

    means = zeros(length(times), length(order_files))
    orders = zeros(length(times), length(order_files))
    std_means = zeros(length(times), length(order_files))

    nlocal_vals = [match(r"(\d+)\.dat$", x).captures[1] for x in order_files]

    # i = 1
    for i in 1:length(order_files)

        raw_data = reinterpret(Float64, read(folder_path * "/" * eta_folders[f] * "/" * exp_files[i]))
        exp_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

        raw_data = reinterpret(Float64, read(folder_path * "/" * eta_folders[f] * "/" * order_files[i]))
        order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

        means[:, i] = mean(exp_data, 2)
        std_means[:, i] = std(exp_data, 2)
        orders[:, i] = mean(order_data, 2)

    end

    # order_p = plot(times, orders, leg = false, xscale = :log10)
    # exp_p   = plot(times, means, leg   = false, xscale = :log10, yscale = :log10)
    # exp_p   = plot(times, means, yerror = std_means, leg   = false, xscale = :log10, yscale = :log10, size = (1024,720))

    order_p = plot(times, orders, lab = nlocal_vals', xscale = :log10)
    exp_p   = plot(times, means, lab = nlocal_vals', xscale = :log10, yscale = :log10)
    exp_p   = plot(times, means, lab = nlocal_vals', xscale = :log10)
    # exp_p   = plot(times, means, yerror = std_means, xscale = :log10, yscale = :log10, size = (1024,720))

    plot(exp_p, order_p, layout = @layout([a b]), size = [1024, 720])
    # plot(exp_p, order_p, layout = @layout[exp_p order_p])

    savefig("$(homedir())/Google\ Drive/proyecto_martin/imagenes/modelo_cvgn/expansion_N_$(N)_$(eta_folders[f]).png")


end
# gui()
