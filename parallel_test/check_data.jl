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

gr()
pyplot()
gui()

### ================================== ###

N = 4096
N = 1024
N = 512
N = 256
N = 128
N = 100
N = 64

τ = 6
τ = 5
τ = 4
τ = 3

times = get_times(0,6)[2:end]

v0 = 1.0

### ================================== ###

folder = "new/NLOC_MET_3D_EXT"
folder = "new/NLOC_P_TOP_3D"

folder = "NLOC_MET_3D_EXT"
folder = "NLOC_P_TOP_3D"

folder_path = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)"
folder_path = "$(homedir())/art_DATA/$(folder)/EXP_N/exp_data_N_$(N)"

eta_folders = readdir(folder_path)

f = 1
data_path = folder_path * "/" * eta_folders[f]

reps = [match(r"\w+(\d+).\w+", x).captures[1]  for x in filter(x -> ismatch(r"^pos_", x), readdir(data_path))]

# r = rand(reps)
r = 1

raw_data = reinterpret(Float64, read(data_path * "/pos_$(r).dat"))

# 3D
pos_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))

x = view(pos_data, :, 1:3:3N)
y = view(pos_data, :, 2:3:3N)
z = view(pos_data, :, 3:3:3N)

# plot(x, y, z, leg = false, size = (800,800), tickfont = font(1), aspect_ratio = :equal)
o = plot(x, y, z, leg = false, size = (800,800))

### ================================== ###

folder = "NLOC_T_MET_3D"
folder = "NLOC_T_TOP_3D"

folder_path = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)"
folder_path = "$(homedir())/art_DATA/$(folder)/EXP/exp_data_N_$(N)"

eta_folders = readdir(folder_path)

f = 1
data_path = folder_path * "/" * eta_folders[f]

reps = [match(r"\w+(\d+).\w+", x).captures[1]  for x in filter(x -> ismatch(r"^pos_", x), readdir(data_path))]
r = 2

raw_data = reinterpret(Float64, read(data_path * "/pos_$(r).dat"))

# 3D
pos_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))

x = view(pos_data, :, 1:3:3N)
y = view(pos_data, :, 2:3:3N)
z = view(pos_data, :, 3:3:3N)

# plot(x, y, z, leg = false, size = (800,800), tickfont = font(1), aspect_ratio = :equal)
e = plot(x, y, z, leg = false, size = (800,800))

### ================================== ###

plot(o, e)

### ================================== ###

size(pos_data)

div(length(raw_data), 3N)

### ================================== ###

order_files = filter( x -> ismatch(r"^order.", x), readdir(folder_path))
exp_files = filter( x -> ismatch(r"^exp.", x), readdir(folder_path))
nn_files = filter( x -> ismatch(r"^nn_mean.", x), readdir(folder_path))

means = zeros(length(times), length(order_files))
orders = zeros(length(times), length(order_files))
nn_means = zeros(length(times), length(order_files))

std_means = zeros(length(times), length(order_files))

### ================================== ###
# para VSK
capt = [match(r"_(\d+\.\d+)\.\w+$|_(\d+\.\d+e-5)\.\w+$", f).captures for f in order_files]
vals = [parse(Float64, vcat(capt...)[i]) for i in find(x -> x != nothing, vcat(capt...))]

# para NLOC
vals = [ parse(Float64, match(r"^\w+_(\d+\.\d+)", x).captures[1]) for x in order_files ]

# para TFLOCK
vals = [ parse(Float64, match(r"(\d\.\d+)\.\w+$", x).captures[1]) for x in order_files ]

### ================================== ###
# i = 3
for i in sortperm(vals)

    println(i)

    raw_data = reinterpret(Float64, read(folder_path * "/" * exp_files[i]))
    exp_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    raw_data = reinterpret(Float64, read(folder_path * "/" * order_files[i]))
    order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    raw_data = reinterpret(Float64, read(folder_path * "/" * nn_files[i]))
    nn_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    means[:, i] = mean(exp_data, 2)
    std_means[:, i] = std(exp_data, 2)
    orders[:, i] = mean(order_data, 2)
    nn_means[:, i] = mean(nn_data, 2)

end

### ================================== ###

gr(size=(1024,720))

plot(times, orders, xscale = :log10, leg = false, xlabel = "t", ylabel = " psi")
png("order_t_m")

plot(times, means, xscale = :log10, yscale = :log10, leg = false, xlabel = "t", ylabel = "rij")
png("rij_t_m")

plot(times, nn_means, xscale = :log10, yscale = :log10, leg = false, xlabel = "t", ylabel = "rnn")
png("rnn_t_m")

plot([vals[i] for i in sortperm(vals)], orders[end, :], leg = false, m = :o, xscale = :log10, xlabel = "k", ylabel = "psi")
png("order_k_m")

[vals[i] for i in sortperm(vals)]
