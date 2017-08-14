using PyPlot,CollectiveDynamics.DataAnalysis
using GR, CollectiveDynamics.DataAnalysis
using Plots, CollectiveDynamics.DataAnalysis


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

times = get_times(τ)

### ================================== ###

folder = "NLOC_DATA"
folder = "NLOC_DATA_3D"
folder = "NLOC_MET_3D"
folder = "NLOC_TOP_3D"
folder = "NLOC_TOP_3D_MEAN"
folder = "TFLOCK_NLOC_DATA"
folder = "TFLOCK_DATA"
folder = "SVM_GRID_FN_3D"
folder = "SVM_GRID_FN_2D"

folder_path = "$(homedir())/art_DATA/$(folder)/EXP/exp_data_N_$(N)"
folder_path = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)"

# eta_folders = readdir(folder_path)

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

# 1,2,3,6
# [order_files[i] for i in [4,5,7,10,13,15]]
# raw_data = reinterpret(Float64, read(folder_path * "/" * order_files[15]))
# order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

### ================================== ###
i = 15
# for i in sortperm(vals)
for i in [4,5,7,10,13,15]

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

writecsv("3D_met_order", hcat(vals[sortperm(vals)], orders[end, sortperm(vals)]))
writecsv("3D_top_order", hcat(vals[sortperm(vals)], orders[end, sortperm(vals)]))
writecsv("3D_vsk_order", hcat(vals[sortperm(vals)], orders[end, sortperm(vals)]))
writecsv("2D_vsk_order", hcat(vals[sortperm(vals)], orders[end, sortperm(vals)]))

plt[:plot](vals[sortperm(vals)], orders[end, sortperm(vals)])
#
# plot(times, order_data)

### ================================== ###

met_vals = readcsv("3D_met_order.csv")
top_vals = readcsv("3D_top_order.csv")
vsk_vals = readcsv("3D_vsk_order.csv")
vsk_vals_2D = readcsv("2D_vsk_order.csv")

### ================================== ###
### PYPLOT
### ================================== ###

plt[:rc]("text", usetex=true)
plt[:rc]("font", family="serif")
plt[:rc]("font", serif="Palatino")
plt[:rc]("font", serif="Times")
plt[:rc]("font", serif="New Century Schoolbook")

### ================================== ###

k_lim = 0.51

y_l = 0.97 #NLOC_DATA
y_h = 1.0 #NLOC_DATA

fs = 14
ls = 12

fig = plt[:figure](num = 1, dpi = 300, facecolor="w", edgecolor="k")
fig[:set_size_inches](2.4, 1.92, forward = true)
fig[:set_size_inches](2, 2, forward = true)

ax = fig[:add_subplot](111)

###==============###==============###==============###

ax[:plot](met_vals[:, 1], met_vals[:, 2], "-o", color = "#7C98AB", ms = 2, lw = 0.5)
ax[:plot](top_vals[:, 1], top_vals[:, 2], "-^", color = "#F6883D", ms = 2, lw = 0.5)

ax[:plot](vcat(vsk_vals, [1.0 0.98])[:, 1], vcat(vsk_vals, [1.0 0.98])[:, 2], "-s", color = "#423b3b", ms = 2, lw = 0.5)

ax[:plot](fill(0.3, 5), collect(linspace(0.,1.,5)), "--", color = "#0e599f", lw = 0.8)

#7C98AB, #F6883D
plt[:xscale]("log")
# plt[:xscale]("linear")

plt[:xlim](0.00085, 1.2)
plt[:xlim](0.00085, 1.8)

plt[:yticks](collect(linspace(0,1,5)))
plt[:yticks](collect(linspace(0,1,3)))

plt[:xticks]([exp10(-3), exp10(-2), exp10(-1), 0.3 , exp10(0)])
plt[:xticklabels](["10^{-3}", "10^{-2}", "10^{-1}", "0.3", "10^{-3}"])

ax[:text](1.0, 0.1, L"\kappa", ha="center", va="center", size=fs)
ax[:text](0.0025, 0.9, L"\Psi(\kappa)", ha="center", va="center", size=fs)

ax[:text](0.2, 0.02, L"\rho_0", ha="center", va="center", size=0.75*fs)
ax[:text](0.3, -0.2, L"\rho_0", ha="center", va="center", size=0.75*fs)

ax[:text](0.85, 0.08, L"\rho", ha="center", va="center", size=fs)
ax[:text](0.0025, 0.9, L"\Psi(\rho)", ha="center", va="center", size=fs)

plt[:tick_params](which = "both", labelsize = ls, direction = "in", pad = 4)
plt[:tight_layout]()

###==============###==============###==============###

# plt[:savefig]("fase_t.eps", dpi = "figure", format = "eps", bbox_inches = "tight" , pad_inches = 0.1)

fig[:savefig]("fase_test1.eps", dpi = 300, format = "eps", bbox_inches = "tight" , pad_inches = 0.1)

fig[:savefig]("order_kappa.eps", dpi = 300, format = "eps", bbox_inches = "tight" , pad_inches = 0.1)

fig[:savefig]("order_vsk.eps", dpi = 300, format = "eps", bbox_inches = "tight" , pad_inches = 0.1)

###==============###==============###==============###

plt[:clf]()

###==============###==============###==============###

vcat(vsk_vals, [1.0 0.98])
