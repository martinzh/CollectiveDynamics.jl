### ============== ### ============== ### ============== ###
## Compute system expansion through
## mean distance between particles
## (export data)
## Martin Zumaya Hernandez
## 21 / 02 / 2017
### ============== ### ============== ### ============== ###

### ================================== ###

using CollectiveDynamics.DataAnalysis

### ================================== ###

N = parse(Int, ARGS[1])

# N = 100
# N = 1024

### ================================== ###

data_folder_path   = "$(homedir())/art_DATA/NLOC_DATA/DATA/data_N_$(N)"
output_folder_path = "$(homedir())/art_DATA/NLOC_DATA/EXP"

folders = readdir(data_folder_path)

params = [match(r"\w+_\d+(_\w+_\d+.\d+)", f).captures[1] for f in folders]
k_vals = [parse(match(r"data_N_\d+_k_(.*)", f).captures[1]) for f in folders]

try
    mkdir(output_folder_path)
catch error
    println("Parent folder already exists")
end

try
    mkdir(output_folder_path * "/exp_data_N_$(N)")
catch error
    println("Output folder already exists")
end

### ================================== ###

### ================================== ###

# f = 1
for f in 1:length(folders)

    psi   = Array{Float64}[]
    means = Array{Float64}[]

    exp_file   = open(output_folder_path * "/exp_data_N_$(N)" * "/exp" * params[f] * ".dat", "w+")
    order_file = open(output_folder_path * "/exp_data_N_$(N)" * "/order" * params[f] * ".dat", "w+")

    ### ================================== ###

    reps = [match(r"\w+_(\d+).\w+", x).captures[1] for x in filter(x -> ismatch(r"pos_\d+.\w+", x), readdir(data_folder_path * "/" * folders[f]))]

    # r = 1
    for r in reps

        println(r)

        raw_data = reinterpret(Float64,read(data_folder_path * "/" * folders[f] * "/pos_$(r).dat"))

        pos_data = reshape(raw_data, 2N, div(length(raw_data), 2N))

        calc_vect_2D_cm(pos_data)

        push!(means, [mean(calc_rij_2D_vect(pos_data[:, i])) for i in 1:size(pos_data,2)])

        raw_data = reinterpret(Float64,read(data_folder_path * "/" * folders[f] * "/vel_$(r).dat"))

        vel_data = reshape(raw_data, 2N, div(length(raw_data), 2N))

        push!(psi, [norm(mean([[vel_data[i, j], vel_data[i+1, j]] for i in 1:2:2N])) for j in 1:size(vel_data, 2)])

    end

    write(exp_file, hcat(means...))
    write(order_file, hcat(psi...))

    close(exp_file)
    close(order_file)

end
### ================================== ###
# using Plots; gr()
#
# τ = 7
#
# times = get_times(τ)
# num_reps = length(filter(x -> ismatch(r"pos_\d+.\w+", x), readdir(data_folder_path * "/" * folders[f])))
#
# raw_data = reinterpret(Float64, read(output_folder_path * "/exp_k_0.0.dat"))
# data = reshape(raw_data, length(times), num_reps)
#
# raw_data = reinterpret(Float64, read(output_folder_path * "/order_k_0.0.dat"))
# order = reshape(raw_data, length(times), num_reps)
#
# gui()
#
# plot(times, mean(data, 2), leg = false, xscale = :log10, yscale = :log10)
# plot(times, mean(order, 2), leg = false, xscale = :log10)
#
# plot(times, order, leg = false, xscale = :log10)
