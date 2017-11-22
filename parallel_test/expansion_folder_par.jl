### ============== ### ============== ### ============== ###
## Compute system expansion through
## mean distance between particles
## (export data)
## Martin Zumaya Hernandez
## 21 / 02 / 2017
### ============== ### ============== ### ============== ###

addprocs(4)

### ================================== ###

using CollectiveDynamics.DataAnalysis
using Plots, CollectiveDynamics.DataAnalysis

### ================================== ###

function calc_Rij_3D(pos::SharedArray, Rij::SharedArray)

    N = size(Rij, 1)

    @parallel for i in 1:3:3N

        for j in (i+3):3:3N

            k = div(i, 3) + 1
            l = div(j, 3) + 1

            Rij[l, k] = norm([pos[i], pos[i+1], pos[i+2]] - [pos[j], pos[j+1], pos[j+2]])
        end
    end

end

### ================================== ###

function calc_Rij_2D(pos, Rij)

    N = size(Rij, 1)

    @parallel for i in 1:2:2N, j in (i+2):2:2N

        k = div(i, 2) + 1
        l = div(j, 2) + 1

        Rij[l, k] = norm([pos[i], pos[i+1]] - [pos[j], pos[j+1]])
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

folder = "NLOC_MET_3D_EXT"
folder = "NLOC_TOP_3D_EXT"

N = 4096
# N = 4000
# N = 1024
# N = 512
# N = 256
# N = 128
# N = 100

k = "9.0"
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

r = 1
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

### ================================== ###
τ = 6
f_times = get_times(τ)

times = union(f_times[1:369],collect(2exp10(5):exp10(5):exp10(6)))

### ================================== ###

sh_pos = pos_data[:,1]


for i in 1:length(Rij)
    Rij[i] = 0.0
end

pos_data
sh_pos

for i in 1:3N
    sh_pos[i] = pos_data[:,1][i]
end

div(length(sh_pos), 3)

calc_Rij_3D(sh_pos, Rij)

mean(Symmetric(Rij, :L))

plot(times, means, xscale = :log10, yscale = :log10, leg = false, m = :o, ms = 1.0)
plot(times, nn_means, xscale = :log10, yscale = :log10, leg = false, m = :o, ms = 1.0)

times = [convert(Int, exp10(i)) for i in 5:6]

for i in 1:(length(times) - 1)

    for t in (times[i]+1):times[i+1]

        if t % times[i] == 0 || t % div(times[i], exp10(1)) == 0
            println(t)
        end
    end

end

Ti = 5
Tf = 6

function get_times(Ti, Tf)

    times = [convert(Int, exp10(i)) for i in Ti:Tf]

    tau = Int64[]

    for i in 1:(length(times) - 1)

        if i > 1

            for t in (times[i]+1):times[i+1]

                if t % times[i] == 0 || t % times[i-1] == 0
                    push!(tau, t)
                end
            end

        else

            for t in (times[i]+1):times[i+1]

                if t % times[i] == 0
                    push!(tau, t)
                end
            end

        end

    end

    return tau
end

tau
