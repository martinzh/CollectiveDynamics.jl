using Plots, CollectiveDynamics.DataAnalysis

### ================================== ###

folder = ARGS[1]
N = parse(Int, ARGS[2])
k = ARGS[3]
w = ARGS[4]
# d = ARGS[5]

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

# data_folder_path   = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)/data_N_$(N)_k_$(k)_w_$(w)"
data_folder_path   = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)/data_N_$(N)_o_$(k)_a_$(w)"

output_folder_path = "$(homedir())/graphs/$(folder)"

make_dir_from_path(output_folder_path)

# params = "_k_$(k)_w_$(w)"
params = "_o_$(k)_a_$(w)"

reps = [match(r"\w+_(\d+).\w+", x).captures[1] for x in filter(x -> ismatch(r"pos_\d+.\w+", x), readdir(data_folder_path))]

trays_plots = Any[]

gr()
# pyplot()

for r in reps

    println(r)

    raw_data = reinterpret(Float64, read(data_folder_path * "/pos_$(r).dat"))
    pos_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))
    # pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

    x = view(pos_data, :, 1:3:3N)
    y = view(pos_data, :, 2:3:3N)
    z = view(pos_data, :, 3:3:3N)

    push!(trays_plots, plot(x, y, z, leg = false, tickfont = font(3)))

end

plot(trays_plots...,  size = (1800,1800))
savefig(output_folder_path * "/reps_trays_N_$(N)$(params).png")
