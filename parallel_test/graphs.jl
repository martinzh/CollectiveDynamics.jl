
### ================================== ###

using Plots

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

times = get_times(0, 6)

times = union(f_times[1:369],collect(2exp10(5):exp10(5):exp10(6)))

folder = "NLOC_MET_3D_EXT"
folder = "NLOC_TOP_3D_EXT"
N = 4096

folder_path = joinpath(homedir(),"art_DATA",folder,"EXP_N","exp_data_N_$(N)")
folder_path = joinpath(homedir(),"art_DATA",folder,"DATA","data_N_$(N)")

folders = readdir(folder_path)

### ================================== ###

order_files = filter( x -> ismatch(r"^order.", x), readdir(folder_path))
exp_files = filter( x -> ismatch(r"^exp.", x), readdir(folder_path))
nn_files = filter( x -> ismatch(r"^nn_mean.", x), readdir(folder_path))

means = zeros(length(times), length(order_files))
orders = zeros(length(times), length(order_files))
nn_means = zeros(length(times), length(order_files))

vals = [ parse(Float64, match(r"^\w+_(\d+\.\d+)", x).captures[1]) for x in order_files ]

for i in sortperm(vals)

    println(i)

    raw_data = reinterpret(Float64, read(folder_path * "/" * exp_files[i]))
    exp_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    raw_data = reinterpret(Float64, read(folder_path * "/" * order_files[i]))
    order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    raw_data = reinterpret(Float64, read(folder_path * "/" * nn_files[i]))
    nn_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    means[:, i] = mean(exp_data, 2)
    orders[:, i] = mean(order_data, 2)
    nn_means[:, i] = mean(nn_data, 2)

end

gr()

# plot(times[2:end], means[2:end, :], xscale = :log10, yscale = :log10, m = :o, ms = 0.8, label = [repr(v) for v in vals], legend = :topleft)
plot(times[2:end], means[2:end, :], xscale = :log10, yscale = :log10, label = [repr(v) for v in vals], legend = :topleft)

plot(times, nn_means, xscale = :log10, yscale = :log10, m = :o, ms = 0.8, label = [repr(v) for v in vals], legend = :topleft)
plot(times, nn_means, xscale = :log10, yscale = :log10, label = [repr(v) for v in vals], legend = :topleft)

plot(times, orders, xscale = :log10, legend = false)

plot(vals, orders[end, :], m = :o, xscale = :log10, legend = false)
xlims!((0.005,5.0))
xlims!((0.0005,0.5))
