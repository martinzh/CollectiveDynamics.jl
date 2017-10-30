using PyPlot,CollectiveDynamics.DataAnalysis
using GR, CollectiveDynamics.DataAnalysis
using Plots, CollectiveDynamics.DataAnalysis


### ================================== ###

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

folder_path = "$(homedir())/art_DATA/$(folder)/EXP_N/exp_data_N_$(N)"
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
# para NLOC, valores de K
vals = [ parse(Float64, match(r"^\w+_(\d+\.\d+)", x).captures[1]) for x in order_files ]

# para NLOC, valores de w
vals = [ parse(Float64, match(r"(\d+\.\d+)\.dat$", x).captures[1]) for x in order_files ]

### ================================== ###

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

writecsv("3D_met_order", hcat(vals[sortperm(vals)], orders[end, sortperm(vals)]))

writecsv("3D_top_order", hcat(vals[sortperm(vals)], orders[end, sortperm(vals)]))

### ================================== ###

met_vals = readcsv("3D_met_order.csv")
top_vals = readcsv("3D_top_order.csv")

### ================================== ###
### PLOTS
### ================================== ###

gr()
gui()

# plot(met_vals[:, 1], met_vals[:, 2], "g-o",  top_vals[:, 1], top_vals[:, 2], "y-^", xlim = (0.0, 1.0))
plot(hcat(met_vals[2:end, 1], top_vals[2:end, 1]),  hcat(met_vals[2:end, 2], top_vals[2:end, 2]), xscale = :log10, xlims = [exp10(-3), 1.0], marker = [:o :^], ms = 6, leg = false, tickfont = font(101, 0), grid = false, dpi = 100, size = (800, 600))

savefig("test_gr.png")

tickfont = font(GR.FONT_TIMES_ROMAN)
### ================================== ###



### ================================== ###
### GR
### ================================== ###

plot(hcat(met_vals[2:end, 1], top_vals[2:end, 1]),  hcat(met_vals[2:end, 2], top_vals[2:end, 2]), xlim = (exp10(-3), 1.0))

plot(met_vals[2:end, 1], met_vals[2:end, 2], "g-o",  top_vals[2:end, 1], top_vals[2:end, 2], "y-^", xlim = (0.0, 0.6))
polyline(met_vals[2:end, 1], met_vals[2:end, 2])
polyline(top_vals[2:end, 1], top_vals[2:end, 2])
setmarkersize(3)
setscale(GR.OPTION_X_LOG)
settextfontprec(101, 0)
updategks()

### ================================== ###
### PYPLOT
### ================================== ###

plt[:rc]("text", usetex=true)
plt[:rc]("font", family="serif")
plt[:rc]("font", serif="Palatino")

### ================================== ###

k_lim = 1.01

y_l = 0.97 #NLOC_DATA
y_h = 1.0 #NLOC_DATA

f = plt[:Figure](figsize = (8, 10), dpi = 300, tight_layout = true)

plt[:plot](met_vals[:, 1], met_vals[:, 2], "b-o",  top_vals[:, 1], top_vals[:, 2], "g-^", ms = 10)

plt[:xscale]("log")
plt[:tick_params](labelsize = 18, direction = "in")
plt[:tick_params](which = "minor", direction = "in")
plt[:xlim](0.0, k_lim)

plt[:tight_layout]()

plt[:savefig]("test_img.eps", dpi = 300, format = "eps", bbox_inches = "tight")
plt[:savefig]("test_img1", dpi = 300, format = "png", bbox_inches = "tight")

savefig

plt[:clf]()
### ================================== ###
# INSET SUBPLOTS

plt[:axes]([.6, .25, .3, .2])
plt[:tick_params](labelsize = 16, direction = "in")
plt[:xlim]((k_lim, sortperm(vals)[end]))
plt[:ylim](y_l, y_h)
plt[:plot](vals[sortperm(vals)], orders[end, sortperm(vals)], "-o", ms = 10)

### ================================== ###

plt[:plot](linspace(0, 2pi, 20), sin.(linspace(0, 2pi, 20)), "-o", ms = 10)
plt[:xlabel](L"\theta", fontsize=20)
plt[:ylabel](L"\sin (\theta)", fontsize=20)
plt[:grid](false)
plt[:show]()

plt[:clf]()
