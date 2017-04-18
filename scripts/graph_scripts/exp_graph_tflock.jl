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
loc_flag = ARGS[3]

# N = 1024
# N = 256

# τ = 7
# τ = 6

times = get_times(τ)

# loc_flag = "nl"

if loc_flag == "l"
    folder_path = "$(homedir())/art_DATA/TFLOCK_DATA/EXP/exp_data_N_$(N)"
elseif loc_flag == "nl"
    folder_path = "$(homedir())/art_DATA/TFLOCK_NLOC_DATA/EXP/exp_data_N_$(N)"
end

eta_folders = readdir(folder_path)

# readdir(folder_path * "/" * eta_folders[1])
#
# n_nl_vals = unique([match(r"(\d+\.\d+)\.dat$", f).captures[1] for f in readdir(folder_path * "/" * eta_folders[1]) ] )

### ================================== ###### ================================== ###

f = 1
for f in 1:length(eta_folders)

    println(eta_folders[f])

    order_files = filter( x -> ismatch(r"^order.", x), readdir(folder_path * "/" * eta_folders[f]))
    exp_files = filter( x -> ismatch(r"^exp.", x), readdir(folder_path * "/" * eta_folders[f]))

    means = zeros(length(times), length(order_files))
    orders = zeros(length(times), length(order_files))
    std_means = zeros(length(times), length(order_files))

    param_vals = [parse(Float64, match(r"(\d+\.\d+)\.dat$", x).captures[1]) for x in order_files]

    # if loc_flag == 1
    #     param_vals = [match(r"(\d+\.\d+)\.dat$", x).captures[1] for x in order_files]
    # elseif loc_flag == 0
    #     param_vals = [match(r"(\d+)\.dat$", x).captures[1] for x in order_files]
    # end

    # i = 1
    for i in 1:length(order_files)

        println(order_files[i])

        raw_data = reinterpret(Float64, read(folder_path * "/" * eta_folders[f] * "/" * exp_files[i]))
        exp_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))
        # exp_data = reshape(raw_data, div(length(raw_data), length(times)), length(times))

        raw_data = reinterpret(Float64, read(folder_path * "/" * eta_folders[f] * "/" * order_files[i]))
        order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))
        # order_data = reshape(raw_data, div(length(raw_data), length(times)), length(times))


        means[:, i]  = mean(exp_data, 2)
        orders[:, i] = mean(order_data, 2)
        std_means[:, i] = std(exp_data, 2)

    end

    order_p = plot(times, orders, lab = param_vals', xscale = :log10)
    # exp_p   = plot(times, means, lab = param_vals', xscale = :log10)
    exp_p   = plot(times, means, lab = param_vals', xscale = :log10, yscale = :log10)

    fase_ord = plot(param_vals, transpose(orders)[:, end], ylims = [0.0, 1.1], marker = :o, leg = false)
    # fase_ord = plot(param_vals, transpose(orders)[:, end], xlims = [0.0, 0.01], ylims = [0.0, 1.1], marker = :o, leg = false)
    # plot!(fase_ord, param_vals[9:end], transpose(orders)[9:end, end], marker = :o, inset_subplots = [(1, bbox(0.5w,0.5h,0.35w,0.35h))], subplot=2)
    # plot!(fase_ord, param_vals[10:end], transpose(orders)[10:end, end], marker = :o, ylims = [0.95, 1.05], inset_subplots = [(1, bbox(0.5w,0.5h,0.45w,0.2h))], subplot=2, leg = false)

    # exp_p   = plot(times, means, yerror = std_means, leg = false, xscale = :log10, yscale = :log10, size = [1024,720])

    # plot(exp_p, order_p, layout = @layout([a b]), size = [1024, 720])
    plot(exp_p, order_p, fase_ord, layout = @layout([a b c]), size = [1420, 720])

    if loc_flag == "l"
        savefig("$(homedir())/Google\ Drive/proyecto_martin/imagenes/modelo_cvgn/expansion_N_$(N)_$(eta_folders[f]).png")
        # savefig("$(homedir())/figures/modelo_cvgn/expansion_N_$(N)_$(eta_folders[f]).png")
    elseif loc_flag == "nl"
        savefig("$(homedir())/Google\ Drive/proyecto_martin/imagenes/modelo_cvgn/expansion_NLOC_N_$(N)_$(eta_folders[f]).png")
    end

end

### ================================== ###### ================================== ###

# gui()
