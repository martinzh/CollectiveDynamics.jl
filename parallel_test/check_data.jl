using Plots, Polynomials, LaTeXStrings
using Plots
using PyPlot
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

τ = 7
τ = 6
τ = 5
τ = 4
τ = 3

times = get_times(0,τ)[2:end]
x_vals = broadcast(x -> log10(x), times)

v0 = 1.0

### ================================== ###

folder = "NLOC_MET_3D_EXT_LST"
folder = "NLOC_TOP_3D_EXT_LST"

folder = "new/NLOC_MET_3D_EXT"
folder = "new/NLOC_TOP_3D_EXT"
folder = "NLOC_TOP_3D_EXT"

folder = "NLOC_MET_3D_EXT"
folder = "NLOC_P_TOP_3D"

folder = "COUZIN_3D"

folder_path = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)"
folder_path = "$(homedir())/art_DATA/$(folder)/EXP_N/exp_data_N_$(N)"

eta_folders = readdir(folder_path)

readdir(joinpath(folder_path,eta_folders[11]))

f = 1
data_path = folder_path * "/" * eta_folders[f]

reps = [match(r"\w+(\d+).\w+", x).captures[1]  for x in filter(x -> ismatch(r"^pos_", x), readdir(data_path))]

# r = rand(reps)
r = 3

raw_data = reinterpret(Float64, read(data_path * "/pos_$(r).dat"))
raw_data = reinterpret(Float64, read(joinpath(folder_path,eta_folders[9],"pos_$(r).dat")))

o = "1.0"
a = "1.0"

raw_data = reinterpret(Float64, read(joinpath(homedir(),"art_DATA","COUZIN_3D","DATA","data_N_128","data_N_$(N)_o_$(o)_a_$(a)","pos_1.dat")));

# 3D
pos_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)));

x = view(pos_data, :, 1:3:3N);
y = view(pos_data, :, 2:3:3N);
z = view(pos_data, :, 3:3:3N);

# plot(x, y, z, leg = false, size = (800,800), tickfont = font(1), aspect_ratio = :equal)
plot(x, y, z, leg = false, size = (800,800))
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

vals = [ (parse(Float64, match(r"^\w+_(\d+\.\d+)_\w+_(\d+\.\d+)", x).captures[1]), parse(Float64, match(r"^\w+_(\d+\.\d+)_\w+_(\d+\.\d+)", x).captures[2]))  for x in order_files ]

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

pyplot()

o = [first(vals[i]) for i in sortperm(vals)]
a = [last(vals[i]) for i in sortperm(vals)]
p = orders[end, :]

gui()

length(orders[end,:])
collect(0.25:0.25:2.0)

scatter3d(o, a, orders[end,:])
surface(o, a, orders[end,:])
xlabel!("")
ylabel!("")

savefig("fase_couzin.svg")


scatter3d(fill(1.0, 8), collect(0.25:0.25:2.0), fill(1.0, 8))
scatter3d!(collect(0.25:0.25:2.0), fill(1.0, 8), fill(1.0, 8))

### ================================== ###
i = 16
raw_data = reinterpret(Float64, read(folder_path * "/" * order_files[i]))
order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

plot(times, mean(order_data, 2), xscale = :log10)
plot(times, mean(order_data, 2), xscale = :log10, yscale = :log10)

plot(times, orders, xscale = :log10)
plot(times, means, xscale = :log10)
plot(times, nn_means, xscale = :log10)

### ================================== ###

r_vals = broadcast(x -> log10(x), means)
r_nn_vals = broadcast(x -> log10(x), nn_means)

x_l = 400
x_l = 320

r_fit_vals = [polyfit(x_vals[x_l:end], r_vals[x_l:end, i], 1) for i in sortperm(vals)]
r_nn_fit_vals = [polyfit(x_vals[x_l:end], r_nn_vals[x_l:end, i], 1) for i in sortperm(vals)]

### ================================== ###

gr(size=(1024,720))
gui()
plot(collect(1:100), sin.(collect(1:100)))

plot(times, orders, xscale = :log10, leg = false, xlabel = "t", ylabel = " psi")
png("order_t_m")

plot(times, means, xscale = :log10, yscale = :log10, leg = false, xlabel = "t", ylabel = "rij")
plot(times, means, xscale = :log10, yscale = :log10, leg = :topleft, xlabel = "t", ylabel = "rij", label = [repr(vals[i]) for i in sortperm(vals)])
png("rij_t_t_1k")
png("rij_t_m_1k")

plot(times, nn_means, xscale = :log10, yscale = :log10, leg = false, xlabel = "t", ylabel = "rnn")
plot(times, nn_means, xscale = :log10, yscale = :log10, leg = :topleft, xlabel = "t", ylabel = "rnn", label = [repr(vals[i]) for i in sortperm(vals)])
png("rnn_t_t_1k")
png("rnn_t_m_1k")

plot([vals[i] for i in sortperm(vals)][1:end], orders[end, 1:end], leg = false, m = :o, xlabel = "k", ylabel = "psi", xlims = [0.0, 2.1])

plot([vals[i] for i in sortperm(vals)][1:end], orders[end, 1:end], leg = false, m = :o, xscale = :log10, xlabel = "k", ylabel = "psi", xlims = [9exp10(-4), 9.1])
plot([vals[i] for i in sortperm(vals)][1:end], orders[end, :][1:end], leg = false, m = :o, xscale = :log10, xlabel = "k", ylabel = "psi", xlims = [9exp10(-4), 1.1])
png("order_k_m_1k")
png("order_k_t_1k")

writecsv("order_top_4k.csv",hcat([vals[i] for i in sortperm(vals)], orders[end,:]))
writecsv("order_met_4k.csv",hcat([vals[i] for i in sortperm(vals)], orders[end,:]))

plot(vals[sortperm(vals)][2:end], [exp10(r_fit_vals[i][0]) for i in 2:length(r_fit_vals)], m = :o, lab = "r_ij", xscale = :log10, yscale = :log10, xlims = (0.8exp10(-2),10), xlabel = "k", ylabel = "D(k)")
plot!(vals[sortperm(vals)][2:end], [exp10(r_nn_fit_vals[i][0]) for i in 2:length(r_fit_vals)], m = :o, lab = "r_nn", xscale = :log10, yscale = :log10, xlims = (0.8exp10(-2),10), xlabel = "k", ylabel = "D(k)")

plot(vals[sortperm(vals)][2:end], [exp10(r_fit_vals[i][0]) for i in 2:length(r_fit_vals)], m =:p, lab = "r_ij", xscale = :log10, yscale = :log10, xlims = (9exp10(-4),10), xlabel = "k", ylabel = "D(k)")
plot!(vals[sortperm(vals)][2:end], [exp10(r_nn_fit_vals[i][0]) for i in 2:length(r_fit_vals)], m =:p, lab = "r_nn", xscale = :log10, yscale = :log10, xlims = (9exp10(-4),1.1), xlabel = "k", ylabel = "D(k)")
png("diff_k_m_1k")
png("diff_k_t_1k")

pwd()
