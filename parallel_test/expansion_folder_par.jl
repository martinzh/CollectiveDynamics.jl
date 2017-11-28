### ============== ### ============== ### ============== ###
## Compute system expansion through
## mean distance between particles
## (export data)
## Martin Zumaya Hernandez
## 21 / 02 / 2017
### ============== ### ============== ### ============== ###

# addprocs(4)

### ================================== ###

# using Plots
@everywhere using DistributedArrays

### ================================== ###

# function calc_Rij_2D(pos, Rij)
#
#     N = size(Rij, 1)
#
#     @parallel for i in 1:2:2N, j in (i+2):2:2N
#
#         k = div(i, 2) + 1
#         l = div(j, 2) + 1
#
#         Rij[l, k] = norm([pos[i], pos[i+1]] - [pos[j], pos[j+1]])
#     end
#
# end

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

@everywhere function calc_Rij_3D(pos, Rij)

    N = size(Rij, 1)

    for i in 1:3:3N

        for j in (i+3):3:3N

            k = div(i, 3) + 1
            l = div(j, 3) + 1

            Rij[l, k] = norm([pos[i], pos[i+1], pos[i+2]] - [pos[j], pos[j+1], pos[j+2]])
        end
    end

end

### ================================== ###

@everywhere function calc_means(pos, vel, mean, nn_mean, vel_mean, N)

    Rij = zeros(N, N)

    for j in localindexes(mean)[1]

        calc_Rij_3D(pos[:,j], Rij)

        mean[j] = mean(Symmetric(Rij, :L))
        nn_mean[j] = mean(sort(Symmetric(Rij, :L), 1)[2,:])
        vel_mean[j] = norm(mean([[vel_data[i, j], vel_data[i+1, j], vel_data[i+2, j]] for i in 1:3:3N]))

    end
end
### ================================== ###

folder = ARGS[1]
N = parse(Int, ARGS[2])
k = ARGS[3]
# w = ARGS[4]


# folder = "NLOC_MET_3D"
# folder = "NLOC_TOP_3D"
# folder = "NLOC_DATA"
# folder = "SVM_GRID_3D"
# folder = "SVM_GRID_FN_2D"

# folder = "new/NLOC_MET_3D_EXT"
# folder = "new/NLOC_TOP_3D_EXT"

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

w = "0.5"

Ti = 0
Tf = 6

times = get_times(Ti, Tf)

### ================================== ###

data_folder_path = joinpath(homedir(),"art_DATA",folder,"DATA","data_N_$(N)","data_N_$(N)_k_$(k)_w_$(w)")

output_folder_path = joinpath(homedir(),"art_DATA",folder,"EXP_N")

make_dir_from_path(output_folder_path)
make_dir_from_path(output_folder_path * "/exp_data_N_$(N)")

params = "_k_$(k)_w_$(w)"
println(params)

### ================================== ###

exp_file        = open(output_folder_path * "/exp_data_N_$(N)" * "/exp" * params * ".dat", "w+")
order_file      = open(output_folder_path * "/exp_data_N_$(N)" * "/order" * params * ".dat", "w+")
nn_mean_file    = open(output_folder_path * "/exp_data_N_$(N)" * "/nn_mean" * params * ".dat", "w+")
# vel_mean_file = open(output_folder_path * "/exp_data_N_$(N)" * "/vel_mean" * params * ".dat", "w+")

### ================================== ###

reps = [match(r"\w+_(\d+).\w+", x).captures[1] for x in filter(x -> ismatch(r"pos_\d+.\w+", x), readdir(data_folder_path))]

# r = 1
for r in reps

    println("rep=",r)

    raw_data = reinterpret(Float64,read(data_folder_path * "/pos_$(r).dat"))

    # pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

    d_pos = distribute(reshape(raw_data, 3N, div(length(raw_data), 3N)), procs = workers(), dist = [1,length(workers())])

    raw_data = reinterpret(Float64,read(data_folder_path * "/vel_$(r).dat"))

    # vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

    d_vel = distribute(reshape(raw_data, 3N, div(length(raw_data), 3N)), procs = workers(), dist = [1,length(workers())])

    d_mean = dzeros(length(times))
    d_nn_mean = dzeros(length(times))
    d_vel_mean = dzeros(length(times))

    @sync for p in workers()
        @async remotecall_wait(calc_means, p, d_pos, d_vel, d_mean, d_nn_mean, d_vel_mean, N)
    end

    # calc_means(d_pos, d_vel, d_mean, d_nn_mean, d_vel_mean, N)

    println("finished calc_means")

    write(exp_file, convert(Array, d_mean))
    write(nn_mean_file, convert(Array, d_nn_mean))
    write(order_file, convert(Array, d_vel_mean))

end

close(d_pos)
close(d_vel)
close(d_mean)
close(d_nn_mean)
close(d_vel_mean)

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
