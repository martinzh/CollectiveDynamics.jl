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

folder = ARGS[1]
N      = parse(Int, ARGS[2])
k      = ARGS[3]
model  = ARGS[4]

# folder = "NLOC_MET_3D"
# folder = "NLOC_TOP_3D"
# folder = "NLOC_DATA"
# folder = "SVM_GRID_3D"
# folder = "SVM_GRID_FN_2D"

# N = 4000
# N = 1024
# N = 512
# N = 256
# N = 128
# N = 100

# k = "9.0"
# k = "0.5"
# k = "0.001"


### ================================== ###

if model == "lnl"
    data_folder_path   = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)/data_N_$(N)_k_$(k)_w_0.5"
    params = "_k_$(k)_w_0.5"
elseif model == "svm"
    data_folder_path   = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)/data_N_$(N)_rho_$(k)"
    params = "_rho_$(k)"
end

output_folder_path = "$(homedir())/art_DATA/$(folder)/EXP"

folders = readdir(data_folder_path)

# params = [match(r"\w+_\d+(_\w+_\d+.\d+)", f).captures[1] for f in folders]
# params = [match(r"\w+_\d+(_\w+_\d+.\d+_\w_\d+\.\d+)", f).captures[1] for f in folders]
# k_vals = [parse(match(r"data_N_\d+_k_(\d+\.\d+)_.", f).captures[1]) for f in folders]

make_dir_from_path(output_folder_path)
make_dir_from_path(output_folder_path * "/exp_data_N_$(N)")

### ================================== ###

### ================================== ###

println(params)

psi = Array{Float64}[]
# vel_means = Array{Float64}[]

order_file = open(output_folder_path * "/exp_data_N_$(N)" * "/order" * params * ".dat", "w+")
# vel_mean_file = open(output_folder_path * "/exp_data_N_$(N)" * "/vel_mean" * params * ".dat", "w+")

### ================================== ###

reps = [match(r"\w+_(\d+).\w+", x).captures[1] for x in filter(x -> ismatch(r"pos_\d+.\w+", x), readdir(data_folder_path))]

Rij = zeros(Float64, N, N)

# r = 1
for r in reps

    println(r)

    raw_data = reinterpret(Float64,read(data_folder_path * "/vel_$(r).dat"))

    vel_data = reshape(raw_data, 2N, div(length(raw_data), 2N))
    push!(psi, [norm(mean([[vel_data[i, j], vel_data[i+1, j]] for i in 1:2:2N])) for j in 1:size(vel_data, 2)])
    # push!(vel_means, vcat([mean([[vel_data[i, j], vel_data[i+1, j]] for i in 1:2:2N]) for j in 1:size(vel_data, 2)]...))

    # vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))
    # push!(psi, [norm(mean([[vel_data[i, j], vel_data[i+1, j], vel_data[i+2, j]] for i in 1:3:3N])) for j in 1:size(vel_data, 2)])
    # push!(vel_means, vcat([mean([[vel_data[i, j], vel_data[i+1, j], vel_data[i+2, j]] for i in 1:3:3N]) for j in 1:size(vel_data, 2)]...))

    # println("pass vel and psi")

end

write(order_file, hcat(psi...))
# write(vel_mean_file, hcat(vel_means...))

close(order_file)
# close(vel_mean_file)

println("Done")

### ================================== ###
