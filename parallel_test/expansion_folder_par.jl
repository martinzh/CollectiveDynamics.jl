### ============== ### ============== ### ============== ###
## Compute system expansion through
## mean distance between particles
## (export data)
## Martin Zumaya Hernandez
## 21 / 02 / 2017
### ============== ### ============== ### ============== ###

# addprocs(4)

### ================================== ###

using Plots, CollectiveDynamics.DataAnalysis

### ================================== ###

function calc_Rij_3D(pos::SharedArray, Rij::SharedArray)

    @parallel for i in 1:3:length(pos)

        ri = div(i,3) + 1

        for j in (i+3):3:length(pos)

            rj = div(j,3) + 1

            Rij[rj,ri] = (pos[i]-pos[j])^2 + (pos[i+1]-pos[j+1])^2 + (pos[i+2]-pos[j+2])^2
        end
    end

end

### ================================== ###

function calc_Rij_2D(pos, Rij)

    @parallel for i in 1:2:length(pos)

        ri = div(i,2) + 1

        for j in (i+2):2:length(pos)

            rj = div(j,2) + 1

            Rij[rj,ri] = (pos[i]-pos[j])^2 + (pos[i+1]-pos[j+1])^2
        end
    end

end
### ================================== ###

folder = ARGS[1]
N = parse(Int, ARGS[2])
k = ARGS[3]
# w = ARGS[4]

w = "0.5"
# d = ARGS[5]

# folder = "NLOC_MET_3D"
# folder = "NLOC_TOP_3D"
# folder = "NLOC_DATA"
# folder = "SVM_GRID_3D"
# folder = "SVM_GRID_FN_2D"

# folder = "NLOC_MET_3D_EXT"
# folder = "NLOC_TOP_3D_EXT"

# N = 4096
# N = 4000
# N = 1024
# N = 512
# N = 256
# N = 128
# N = 100

# k = "9.0"
# k = "0.5"
# k = "0.001"
# k = "0.0125"
# w = "0.5"

# τ = 6
# times = get_times(τ)

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

### ================================== ###

println(params)

# psi       = Array{Float64}[]
# means     = Array{Float64}[]
# nn_means  = Array{Float64}[]
# vel_means = Array{Float64}[]

exp_file   = open(output_folder_path * "/exp_data_N_$(N)" * "/exp" * params * ".dat", "w+")
nn_mean_file = open(output_folder_path * "/exp_data_N_$(N)" * "/nn_mean" * params * ".dat", "w+")

order_file = open(output_folder_path * "/exp_data_N_$(N)" * "/order" * params * ".dat", "w+")
# vel_mean_file = open(output_folder_path * "/exp_data_N_$(N)" * "/vel_mean" * params * ".dat", "w+")

### ================================== ###

reps = [match(r"\w+_(\d+).\w+", x).captures[1] for x in filter(x -> ismatch(r"pos_\d+.\w+", x), readdir(data_folder_path))]

nn_means = []
means    = []
psi      = Array{Float64}[]

Rij = SharedArray{Float64}(N,N)
sh_pos = SharedArray{Float64}(3N)

# r = 3
for r in reps

    println(r)

    raw_data = reinterpret(Float64,read(data_folder_path * "/pos_$(r).dat"))
    pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

    raw_data = reinterpret(Float64,read(data_folder_path * "/vel_$(r).dat"))
    vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

    for j in 1:size(pos_data,2)

        println(j)

        for i in 1:3N
            sh_pos[i] = pos_data[:,j][i]
        end

        calc_Rij_3D(sh_pos, Rij)

        push!(means,mean(Symmetric(Rij, :L)))
        push!(nn_means, mean([sort(Symmetric(fetch(Rij), :L)[:, k])[2] for k in 1:size(Rij,2)]))

    end

    push!(psi, [norm(mean([[vel_data[i, j], vel_data[i+1, j], vel_data[i+2, j]] for i in 1:3:3N])) for j in 1:size(vel_data, 2)])

    println("pass vel and psi")

end


# write(exp_file, hcat(means...))
# write(nn_mean_file, hcat(nn_means...))

write(exp_file, means)
write(nn_mean_file, nn_means)

write(order_file, hcat(psi...))
# write(vel_mean_file, hcat(vel_means...))

close(exp_file)
close(nn_mean_file)

close(order_file)
# close(vel_mean_file)

println("Done")

### ================================== ###

# gui()
# plot(union(times[1:369],collect(2exp10(5):exp10(5):exp10(6))), means, m = :o, xscale = :log10, yscale = :log10, leg = false)
# plot(union(times[1:369],collect(2exp10(5):exp10(5):exp10(6))), nn_means, m = :o, xscale = :log10, yscale = :log10, leg = false)
#
# plot(union(times[1:369],collect(2exp10(5):exp10(5):exp10(6))), [norm(mean([[vel_data[i, j], vel_data[i+1, j], vel_data[i+2, j]] for i in 1:3:3N])) for j in 1:size(vel_data, 2)], xscale = :log10, leg = false)
