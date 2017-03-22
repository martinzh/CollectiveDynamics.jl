### ============== ### ============== ### ============== ###
## Create and export order and expansion graphs
## (inertial model)
## Martin Zumaya Hernandez
## 22 / 03 / 2017
### ============== ### ============== ### ============== ###

### ================================== ###### ================================== ###

using Plots, CollectiveDynamics.DataAnalysis

# gui()

# gr()
pyplot()

### ================================== ###### ================================== ###

N = parse(Int, ARGS[1])
τ = parse(Int, ARGS[2])
loc_flag = parse(Int, ARGS[3])

# N = 1024
# N = 256

# τ = 7
# τ = 6

times = get_times(τ)

loc_flag = 1

if loc_flag == 1
    folder_path = "$(homedir())/art_DATA/TFLOCK_DATA/EXP/exp_data_N_$(N)"
elseif loc_flag == 0
    folder_path = "$(homedir())/art_DATA/TFLOCK_NLOC_DATA/EXP/exp_data_N_$(N)"
end

eta_folders = readdir(folder_path)

### ================================== ###### ================================== ###

# f = 1
for f in 1:length(eta_folders)

    println(eta_folders[f])

    order_files = filter( x -> ismatch(r"^order.", x), readdir(folder_path * "/" * eta_folders[f]))
    exp_files = filter( x -> ismatch(r"^exp.", x), readdir(folder_path * "/" * eta_folders[f]))

    means = zeros(length(times), length(order_files))
    orders = zeros(length(times), length(order_files))
    std_means = zeros(length(times), length(order_files))

    if loc_flag == 1
        param_vals = [match(r"(\d+\.\d+)\.dat$", x).captures[1] for x in order_files]
    elseif loc_flag == 0
        paam_vals = [match(r"(\d+)\.dat$", x).captures[1] for x in order_files]
    end

    # i = 1
    for i in 1:length(order_files)

        raw_data = reinterpret(Float64, read(folder_path * "/" * eta_folders[f] * "/" * exp_files[i]))
        exp_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

        raw_data = reinterpret(Float64, read(folder_path * "/" * eta_folders[f] * "/" * order_files[i]))
        order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

        means[:, i]  = mean(exp_data, 2)
        orders[:, i] = mean(order_data, 2)
        std_means[:, i] = std(exp_data, 2)

    end

    order_p = plot(times, orders, lab = param_vals', xscale = :log10)
    exp_p   = plot(times, means, lab = param_vals', xscale = :log10, yscale = :log10)

    # exp_p   = plot(times, means, yerror = std_means, leg   =© false, xscale = :log10, yscale = :log10, size = [1024,720])

    plot(exp_p, order_p, layout = @layout([a b]), size = [1024, 720])

    savefig("$(homedir())/Google\ Drive/proyecto_martin/imagenes/modelo_cvgn/expansion_N_$(N)_$(eta_folders[f]).png")

end

### ================================== ###### ================================== ###
