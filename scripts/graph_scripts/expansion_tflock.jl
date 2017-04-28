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
function calc_Rij(pos, Rij)

    N = size(Rij, 1)

    for i in 1:3:3N, j in (i+3):3:3N

        k = div(i, 3) + 1
        l = div(j, 3) + 1

        Rij[k, l] = norm([pos[i], pos[i+1], pos[i+2]] - [pos[j], pos[j+1], pos[j+2]])
        Rij[l, k] = Rij[k, l]
    end

end
### ================================== ###

N = parse(Int, ARGS[1])
loc_flag = parse(Int, ARGS[2])

# N = 100
# N = 1024
# N = 128
# N = 256

# loc_flag = 1

v0 = 0.1
### ================================== ###

if loc_flag == 1
    data_folder_path   = "$(homedir())/art_DATA/TFLOCK_DATA/DATA/data_N_$(N)"
    output_folder_path = "$(homedir())/art_DATA/TFLOCK_DATA/EXP"
elseif loc_flag == 0
    data_folder_path   = "$(homedir())/art_DATA/TFLOCK_NLOC_DATA/DATA/data_N_$(N)"
    output_folder_path = "$(homedir())/art_DATA/TFLOCK_NLOC_DATA/EXP"
end

make_dir_from_path(output_folder_path)

eta_folders = readdir(data_folder_path)
eta_vals = unique([match(r"\w+_(\d+\.\d+)", f).captures[1] for f in eta_folders])

### ================================== ###

# f = 1
for f in 1:length(eta_folders)

    println(eta_folders[f])

    make_dir_from_path(output_folder_path * "/exp_data_N_$(N)")
    make_dir_from_path(output_folder_path * "/exp_data_N_$(N)/eta_$(eta_vals[f])")

    noise_folders = readdir(data_folder_path * "/" * eta_folders[f])

    # i = 1
    for i in 1:length(noise_folders)

        println(noise_folders[i])

        psi   = Array{Float64}[]
        means = Array{Float64}[]
        nn_means = Array{Float64}[]

        exp_file   = open(output_folder_path * "/exp_data_N_$(N)/eta_$(eta_vals[f])" * "/exp_" * noise_folders[i] * ".dat", "w+")
        order_file = open(output_folder_path * "/exp_data_N_$(N)/eta_$(eta_vals[f])" * "/order_" * noise_folders[i] * ".dat", "w+")
        nn_mean_file = open(output_folder_path * "/exp_data_N_$(N)" * "/nn_mean" * params[f] * ".dat", "w+")

        ### ================================== ###

        reps = [match(r"\w+_(\d+).\w+", x).captures[1] for x in filter(x -> ismatch(r"pos_\d+.\w+", x), readdir(data_folder_path * "/" * eta_folders[f] * "/" * noise_folders[i]))]

        Rij = zeros(Float64, N, N)

        # r = 1
        for r in reps

            println(r)

            raw_data = reinterpret(Float64,read(data_folder_path * "/" * eta_folders[f] * "/" * noise_folders[i] * "/" * "/pos_$(r).dat"))
            pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

            calc_vect_3D_cm(pos_data)

            push!(means, [mean(calc_rij_3D_vect(pos_data[:, i])) for i in 1:size(pos_data,2)])

            for j in 1:size(pos_data, 2)

                calc_Rij(pos_data[:, j], Rij)
                nn_time[j] = mean([sort(Rij[:, i])[2] for i in 1:N])

            end

            push!(nn_means, nn_time)

            raw_data = reinterpret(Float64,read(data_folder_path * "/" * eta_folders[f] * "/" * noise_folders[i] * "/" * "/vel_$(r).dat"))
            vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

            push!(psi, (1. / v0) * [norm(mean([[vel_data[i, j], vel_data[i+1, j], vel_data[i+2, j]] for i in 1:3:3N])) for j in 1:size(vel_data, 2)])

        end

        write(exp_file, hcat(means...))
        write(order_file, hcat(psi...))
        write(nn_mean_file, hcat(nn_means...))

        close(exp_file)
        close(order_file)
        close(nn_mean_file)

    end

end

### ================================== ###
