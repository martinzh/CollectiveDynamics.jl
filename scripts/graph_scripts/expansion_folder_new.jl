### ============== ### ============== ### ============== ###
## Compute system expansion through
## mean distance between particles
## (export data)
## Martin Zumaya Hernandez
## 21 / 02 / 2017
### ============== ### ============== ### ============== ###

### ================================== ###

# using CollectiveDynamics.DataAnalysis

### ================================== ###

function make_dir_from_path(path)

    try
        mkdir(path)
    catch error
        println("Main data folder already exists")
    end

end

### ================================== ###

function get_times(Ti, Tf)

    times = [convert(Int, exp10(i)) for i in Ti:Tf]
    tau = Int64[]

    push!(tau, 1)

    for i in 1:(length(times) - 1)

        for t in (times[i]+1):times[i+1]

            if t % times[i] == 0 || t % div(times[i], exp10(1)) == 0
                push!(tau, t)
                # println("//////// ", t)
            end
        end

    end

    return tau
end

### ================================== ###
function calc_Rij_3D(pos, Rij)

    N = size(Rij, 1)

    for i in 1:3:3N, j in (i+3):3:3N

        k = div(i, 3) + 1
        l = div(j, 3) + 1

        Rij[k, l] = norm([pos[i], pos[i+1], pos[i+2]] - [pos[j], pos[j+1], pos[j+2]])
        # Rij[l, k] = Rij[k, l]
    end

end
### ================================== ###
function calc_Rij_2D(pos, Rij)

    N = size(Rij, 1)

    for i in 1:2:2N, j in (i+2):2:2N

        k = div(i, 2) + 1
        l = div(j, 2) + 1

        Rij[k, l] = norm([pos[i], pos[i+1]] - [pos[j], pos[j+1]])
        # Rij[l, k] = Rij[k, l]
    end

end
### ================================== ###

folder = ARGS[1]
N = parse(Int, ARGS[2])
k = ARGS[3]
w = ARGS[4]
d = ARGS[5]

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

data_folder_path   = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)/data_N_$(N)_k_$(k)_w_$(w)"
# data_folder_path   = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)/eta_1.5/eta_1.5_T_0.01_nl_$(k)"
# data_folder_path   = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)/data_N_$(N)_rho_$(k)"

output_folder_path = "$(homedir())/art_DATA/$(folder)/EXP_N"

folders = readdir(data_folder_path)

params = "_k_$(k)_w_$(w)"
# params = "_eta_1.5_T_0.01_n_l_$(k)"
# params = "_rho_$(k)"

# params = [match(r"\w+_\d+(_\w+_\d+.\d+)", f).captures[1] for f in folders]
# params = [match(r"\w+_\d+(_\w+_\d+.\d+_\w_\d+\.\d+)", f).captures[1] for f in folders]
# k_vals = [parse(match(r"data_N_\d+_k_(\d+\.\d+)_.", f).captures[1]) for f in folders]

make_dir_from_path(output_folder_path)
make_dir_from_path(output_folder_path * "/exp_data_N_$(N)")

### ================================== ###

times = get_times(0,6)

### ================================== ###

println(params)

psi       = Array{Float64}[]
means     = Array{Float64}[]
nn_means  = Array{Float64}[]
# vel_means = Array{Float64}[]

exp_file   = open(output_folder_path * "/exp_data_N_$(N)" * "/exp" * params * ".dat", "w+")
nn_mean_file = open(output_folder_path * "/exp_data_N_$(N)" * "/nn_mean" * params * ".dat", "w+")

order_file = open(output_folder_path * "/exp_data_N_$(N)" * "/order" * params * ".dat", "w+")
# vel_mean_file = open(output_folder_path * "/exp_data_N_$(N)" * "/vel_mean" * params * ".dat", "w+")

### ================================== ###

reps = [match(r"\w+_(\d+).\w+", x).captures[1] for x in filter(x -> ismatch(r"pos_\d+.\w+", x), readdir(data_folder_path))]

Rij = zeros(Float64, N, N)

# r = 1
for r in reps

    println(r)

    raw_data = reinterpret(Float64,read(data_folder_path * "/pos_$(r).dat"))

    # if d == "2"
    #     pos_data = reshape(raw_data, 2N, div(length(raw_data), 2N))
    #     # calc_vect_2D_cm(pos_data)
    #     push!(means, [mean(calc_rij_2D_vect(pos_data[:, i])) for i in 1:size(pos_data,2)])
    # elseif d == "3"
    #     pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))
    #     # calc_vect_3D_cm(pos_data)
    #     push!(means, [mean(calc_rij_3D_vect(pos_data[:, i])) for i in 1:size(pos_data,2)])
    # end

    if d == "2"
        pos_data = reshape(raw_data, 2N, div(length(raw_data), 2N))
        # calc_vect_2D_cm(pos_data)
        push!(means, [mean(calc_rij_2D_vect(pos_data[:, i])) for i in 1:size(pos_data,2)])

    elseif d == "3"

        if length(raw_data) == 3N*(length(times)-1)
            pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))
        else
            pos_data = reshape(raw_data[3N+1:end], 3N, div(length(raw_data[3N+1:end]), 3N))
        end

        # pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))
        # calc_vect_3D_cm(pos_data)

        push!(means, [mean(calc_rij_3D_vect(pos_data[:, i])) for i in 1:size(pos_data,2)])
    end

    # println("pass pos_data")

    nn_time = zeros(size(pos_data, 2))

    for j in 1:size(pos_data, 2)

        if d == "2" calc_Rij_2D(pos_data[:, j], Rij) end
        if d == "3" calc_Rij_3D(pos_data[:, j], Rij) end

        # nn_time[j] = mean([sort(Rij[:, i])[2] for i in 1:N])

        nn_mean[j] = mean(sort(Symmetric(Rij, :L), 1)[2,:])

    end

    push!(nn_means, nn_time)

    println("pass nn_means")

    raw_data = reinterpret(Float64,read(data_folder_path * "/vel_$(r).dat"))

    # if d == "2"
    #     vel_data = reshape(raw_data, 2N, div(length(raw_data), 2N))
    #
    #     push!(psi, [norm(mean([[vel_data[i, j], vel_data[i+1, j]] for i in 1:2:2N])) for j in 1:size(vel_data, 2)])
    #     # push!(vel_means, vcat([mean([[vel_data[i, j], vel_data[i+1, j]] for i in 1:2:2N]) for j in 1:size(vel_data, 2)]...))
    #
    # elseif d == "3"
    #     vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))
    #
    #     push!(psi, [norm(mean([[vel_data[i, j], vel_data[i+1, j], vel_data[i+2, j]] for i in 1:3:3N])) for j in 1:size(vel_data, 2)])
    #     # push!(vel_means, vcat([mean([[vel_data[i, j], vel_data[i+1, j], vel_data[i+2, j]] for i in 1:3:3N]) for j in 1:size(vel_data, 2)]...))
    # end

    if d == "2"
        vel_data = reshape(raw_data, 2N, div(length(raw_data), 2N))

        push!(psi, [norm(mean([[vel_data[i, j], vel_data[i+1, j]] for i in 1:2:2N])) for j in 1:size(vel_data, 2)])
        # push!(vel_means, vcat([mean([[vel_data[i, j], vel_data[i+1, j]] for i in 1:2:2N]) for j in 1:size(vel_data, 2)]...))

    elseif d == "3"

        if length(raw_data) == 3N*(length(times)-1)
            vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))
        else
            vel_data = reshape(raw_data[3N+1:end], 3N, div(length(raw_data[3N+1:end]), 3N))
        end

        # vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

        push!(psi, [norm(mean([[vel_data[i, j], vel_data[i+1, j], vel_data[i+2, j]] for i in 1:3:3N])) for j in 1:size(vel_data, 2)])
        # push!(vel_means, vcat([mean([[vel_data[i, j], vel_data[i+1, j], vel_data[i+2, j]] for i in 1:3:3N]) for j in 1:size(vel_data, 2)]...))
    end

    println("pass vel and psi")

end

write(exp_file, hcat(means...))
write(nn_mean_file, hcat(nn_means...))

write(order_file, hcat(psi...))
# write(vel_mean_file, hcat(vel_means...))

close(exp_file)
close(nn_mean_file)

close(order_file)
# close(vel_mean_file)

println("Done")

### ================================== ###
