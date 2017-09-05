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

folder = "NLOC_MET_2D"
folder = "NLOC_MET_3D"

folder = "NLOC_TOP_2D"
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

# para TFLOCK
vals = [ parse(Float64, match(r"(\d\.\d+)\.\w+$", x).captures[1]) for x in order_files ]

# 1,2,3,6
# [order_files[i] for i in [4,5,7,10,13,15]]
# raw_data = reinterpret(Float64, read(folder_path * "/" * order_files[15]))
# order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

### ================================== ###
# i = 1
for i in sortperm(vals)

    println(i)

    # raw_data = reinterpret(Float64, read(folder_path * "/" * exp_files[i]))
    # exp_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    raw_data = reinterpret(Float64, read(folder_path * "/" * order_files[i]))
    order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    # raw_data = reinterpret(Float64, read(folder_path * "/" * nn_files[i]))
    # nn_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    orders[:, i] = mean(order_data, 2)
    # means[:, i] = mean(exp_data, 2)
    # std_means[:, i] = std(exp_data, 2)
    # nn_means[:, i] = mean(nn_data, 2)

end

insert!(vals[sortperm(vals)], 8, 0.6)
insert!(orders[end, sortperm(vals)], 8, 0.46)

psi_vals = vcat(hcat(vals[sortperm(vals)], orders[end, sortperm(vals)]), [10.0 0.95; 20.0 0.95])
psi_vals[1,1] = 0.001

writecsv("2D_met_order.csv", psi_vals)
writecsv("2D_top_order.csv", psi_vals)

writecsv("3D_met_order.csv", hcat(vals[sortperm(vals)], orders[end, sortperm(vals)]))
writecsv("3D_top_order.csv", hcat(vals[sortperm(vals)], orders[end, sortperm(vals)]))

writecsv("3D_vsk_order.csv", hcat(vals[sortperm(vals)], orders[end, sortperm(vals)]))
writecsv("2D_vsk_order.csv", hcat(vals[sortperm(vals)], orders[end, sortperm(vals)]))

writecsv("tflock_order.csv", hcat(vals[sortperm(vals)], orders[end, sortperm(vals)]))

### ================================== ###
## TEST

f = plot(vals[sortperm(vals)], orders[end, sortperm(vals)], marker = :o, leg = false, xscale = :log10, xlims = [0.1, 5.5])
plot!(f, fill(0.3, 10), collect(0.0:0.1:1.0))

plot(psi_vals[:,1], psi_vals[:,2], marker = :o, leg = false, xscale = :log10, xlims = [0.0005, 20.5])

plot(vals[sortperm(vals)], orders[end, sortperm(vals)], marker = :o, leg = false, xscale = :log10, xlims = [0.0001, 1.5])
plot(vals[sortperm(vals)], orders[end, sortperm(vals)], leg = false, xscale = :log10, xlims = [0.0001, 1.5])

plot(vals[sortperm(vals)], orders[end, sortperm(vals)], marker = :o, leg = false, xscale = :log10, xlims = [0.01, 1.0])

plot(vals[sortperm(vals)], orders[end, sortperm(vals)], leg = false, xscale = :log10, xlims = [0.01, 1.5])

plot(times, orders, xscale = :log10, leg = false)

#TFLOCK_NLOC_DATA
plot(vals[sortperm(vals)], inv(0.1) * orders[end, sortperm(vals)], marker = :o, leg = false, xscale = :log10, xlims = [0.0001, 1.0])
plot(times, inv(0.1) * orders, xscale = :log10, leg = false)

plot(times, orders, xscale = :lin, leg = false)
gr()

ax[:plot](vals[sortperm(vals)], inv(0.1) .* orders[end, sortperm(vals)], "-o", color = "#000000", ms = 2, lw = 0.5)

### ================================== ###

met_vals_2D = readcsv("2D_met_order.csv")
top_vals_2D = readcsv("2D_top_order.csv")

met_vals = readcsv("3D_met_order.csv")
top_vals = readcsv("3D_top_order.csv")

tflock_vals = readcsv("tflock_order.csv")

vsk_vals = readcsv("3D_vsk_order.csv")
vsk_vals_2D = readcsv("2D_vsk_order.csv")
vsk_vals_2D = readcsv("$(homedir())/art_DATA/order_vicsek_n10000.csv")

### ================================== ###
### PYPLOT
### ================================== ###

plt[:rc]("text", usetex=true)
plt[:rc]("font", family="serif")
plt[:rc]("font", serif="New Century Schoolbook")
# plt[:rc]("font", serif="Palatino")
# plt[:rc]("font", serif="Times")

### ================================== ###

k_lim = 0.51

y_l = 0.97 #NLOC_DATA
y_h = 1.0 #NLOC_DATA

fs = 12
ls = 10

fig = plt[:figure](num = 1, dpi = 300, facecolor="w", edgecolor="k")
fig[:set_size_inches](2.4, 1.92, forward = true)
fig[:set_size_inches](2, 2, forward = true)

ax = fig[:add_subplot](111)

###==============###==============###==============###

plt[:clf]()

###==============###==============###==============###
ax = fig[:add_subplot](111)

ax[:plot](met_vals[:, 1], met_vals[:, 2], "-o", color = "#7C98AB", ms = 2, lw = 0.5, label = "metric")
ax[:plot](top_vals[:, 1], top_vals[:, 2], "-^", color = "#F6883D", ms = 2, lw = 0.5, label = "topological")

ax[:plot](met_vals_2D[:, 1], met_vals_2D[:, 2], "-o", color = "#7C98AB", ms = 2, lw = 0.5, label = "metric")
ax[:plot](vcat(top_vals_2D[:, 1], 10.0), vcat(top_vals_2D[:, 2], 0.953), "-^", color = "#F6883D", ms = 2, lw = 0.5, label = "topological")
ax[:plot](met_vals[:, 1], met_vals[:, 2], "-o", color = "#000000", ms = 2, lw = 0.5)

ax[:plot](tflock_vals[:, 1], inv(0.1) .* tflock_vals[:, 2], "-o", color = "#000000", ms = 2, lw = 0.5)
###==============###==============###==============###

ax[:plot](vcat(vsk_vals, [1.0 0.98])[:, 1], vcat(vsk_vals, [1.0 0.98])[:, 2], "-s", color = "#423b3b", ms = 2, lw = 0.5)

###==============###==============###==============###

ax[:plot](vsk_vals_2D[:, 1], vsk_vals_2D[:, end], "-s", color = "#423b3b", ms = 2, lw = 0.5)

###==============###==============###==============###

# ρ_0
ax[:plot](fill(0.3, 5), collect(linspace(0.,1.,5)), "--", color = "#0e599f", lw = 0.8)

# unbounded
ax[:plot](vcat(vsk_vals, [1.0 0.98])[:, 1], [vsk_vals[1, 2] + rand(-0.005:0.001:0.005) for i in 1:16], ":", color = "#CD2B18", lw = 1.2)

# unbounded
ax[:plot](vsk_vals_2D[:, 1], [vsk_vals_2D[1, end] + rand(-0.005:0.001:0.005) for i in 1:size(vsk_vals_2D, 1)], ":", color = "#CD2B18", lw = 1.2)

#7C98AB, #F6883D
plt[:xscale]("log")
# plt[:xscale]("linear")

met_vals[1,1] = 0.001
top_vals[1,1] = 0.0001

println(met_vals)
println(top_vals)

plt[:xlim](0.0005, 15)
plt[:xlim](0.00085, 1.8)
plt[:xlim](0.00008, 0.51)
plt[:xlim](0.00005, 30)
plt[:xlim](0.000008, 0.6)

plt[:ylim](0., 1.01)
plt[:ylim](0.2, 1.05) # tflock

plt[:yticks](collect(linspace(0,1,5)))
plt[:yticks](collect(linspace(0,1,3)))
plt[:yticks](collect(linspace(0.2,1,3))) # tflock

plt[:xticks]([exp10(-4), exp10(-2), exp10(0)])
plt[:xticks]([exp10(-4), exp10(-2), exp10(0), exp10(1)])

plt[:xticks]([exp10(-3), exp10(-2), exp10(-1), 0.3 , exp10(0)])
plt[:xticklabels](["10^{-3}", "10^{-2}", "10^{-1}", "0.3", "10^{-3}"])

plt[:xticks]([exp10(-4), exp10(-3), exp10(-2), exp10(-1)])

# 2D MET, TOP
ax[:text](18, 0.05, L"\kappa", ha="center", va="center", size=fs)
ax[:text](2exp10(-4), 0.9, L"\Psi(\kappa)", ha="center", va="center", size=fs)

# TFLOCK
ax[:text](0.4, 0.3, L"\kappa", ha="center", va="center", size=fs)
ax[:text](2.8exp10(-5), 0.97, L"\Psi(\kappa)", ha="center", va="center", size=fs)

# 3D
ax[:text](11, 0.05, L"\kappa", ha="center", va="center", size=fs)
ax[:text](0.00025, 0.94, L"\Psi(\kappa)", ha="center", va="center", size=fs)

ax[:text](0.3, -0.2, L"\rho_0", ha="center", va="center", size=0.75*fs) # 3D
ax[:text](0.2, 0.08, L"\rho_0", ha="center", va="center", size=0.75*fs) # 2D

ax[:text](0.85, 0.08, L"\rho", ha="center", va="center", size=fs) # 3D
ax[:text](0.6, -0.12, L"\rho", ha="center", va="center", size=fs) # 2D

ax[:text](0.0025, 0.9, L"\Psi(\rho)", ha="center", va="center", size=fs) # 3D
ax[:text](0.00025, 0.95, L"\Psi(\rho)", ha="center", va="center", size=fs) # 2D
ax[:tick_params](axis="x",which="minor",bottom="off")
ax[:tick_params](axis="y",which="minor",left="off")

plt[:tick_params](which = "both", labelsize = ls, direction = "in", pad = 4)
plt[:tight_layout]()

ax[:set_xlabel](L"\kappa", labelpad =0)
ax[:set_ylabel](L"\Psi", labelpad =0)

plt[:legend](bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0., fontsize = 8.5)

###==============###==============###==============###

# plt[:savefig]("fase_t.eps", dpi = "figure", format = "eps", bbox_inches = "tight" , pad_inches = 0.1)

fig[:savefig]("fase_test1.eps", dpi = 300, format = "eps", bbox_inches = "tight" , pad_inches = 0.1)

fig[:savefig]("order_kappa.eps", dpi = 300, format = "eps", bbox_inches = "tight" , pad_inches = 0.1)
fig[:savefig]("order_kappa.png", dpi = 300, format = "png", bbox_inches = "tight" , pad_inches = 0.1)

fig[:savefig]("order_vsk.eps", dpi = 300, format = "eps", bbox_inches = "tight" , pad_inches = 0.1)

fig[:savefig]("2D_order_vsk.eps", dpi = 300, format = "eps", bbox_inches = "tight" , pad_inches = 0.1)

fig[:savefig]("met_order_kappa.png", dpi = 300, format = "png", bbox_inches = "tight" , pad_inches = 0.1)

fig[:savefig]("2D_order_kappa.eps", dpi = 300, format = "eps", bbox_inches = "tight" , pad_inches = 0.1)
fig[:savefig]("2D_order_kappa.png", dpi = 300, format = "png", bbox_inches = "tight" , pad_inches = 0.1)

fig[:savefig]("tflock_order_kappa.eps", dpi = 300, format = "eps", bbox_inches = "tight" , pad_inches = 0.1)
fig[:savefig]("tflock_order_kappa.png", dpi = 300, format = "png", bbox_inches = "tight" , pad_inches = 0.1)

###==============###==============###==============###

plt[:clf]()

###==============###==============###==============###
