ENV["PLOTS_USE_ATOM_PLOTPLANE"] = "false"

using PyPlot, LaTeXStrings, CollectiveDynamics.DataAnalysis

###==============###==============###==============###

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

tau = get_times(τ)

v0 = 1.0
###==============###==============###==============###

folder = "NLOC_DATA"
folder = "NLOC_DATA_3D"

folder = "NLOC_MET_2D"
folder = "NLOC_TOP_2D"

folder = "NLOC_MET_3D"
folder = "NLOC_TOP_3D"

folder = "NLOC_TOP_3D_MEAN"

folder = "TFLOCK_NLOC_DATA"
folder = "TFLOCK_DATA"

folder_path = "$(homedir())/art_DATA/$(folder)/EXP/exp_data_N_$(N)"
folder_path = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)"
folder_path = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)/eta_1.5"

eta_folders = readdir(folder_path)

f = 3
data_path = folder_path * "/" * eta_folders[f]

reps = [match(r"\w+(\d+).\w+", x).captures[1]  for x in filter(x -> ismatch(r"^pos_", x), readdir(data_path))]

r = rand(reps)
r = 8

raw_data = reinterpret(Float64, read(data_path * "/pos_$(r).dat"))

###==============###==============###==============###

pos_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))

x = view(pos_data, :, 1:3:3N)
y = view(pos_data, :, 2:3:3N)
z = view(pos_data, :, 3:3:3N)

###==============###==============###==============###

pos_data = transpose(reshape(raw_data, 2N, div(length(raw_data), 2N)))

x = view(pos_data, :, 1:2:2N)
y = view(pos_data, :, 2:2:2N)


###==============###==============###==============###

plt[:rc]("text", usetex=true)
plt[:rc]("font", family="serif")
plt[:rc]("font", serif="New Century Schoolbook")

# plt[:rc]("figure.autolayout")
# plt[:rcParams]["figure.autolayout"] = true

plt[:clf]()

###==============###==============###==============###
fs = 8
ls = 6

fs = 10
ls = 8

sx = sy = 6

###==============###==============###==============###

plt[:ticklabel_format](style = "sci", scilimits = (2,4))
plt[:ticklabel_format](style = "plain")

fig = plt[:figure](num = 1, dpi = 100, facecolor="w", edgecolor="k")
fig[:set_size_inches](sx, sy, forward = true)

fig = plt[:figure](num = 1, dpi = 300, facecolor="w", edgecolor="w")
fig[:set_size_inches](2.4, 1.92, forward = true)

# fig[:subplots_adjust](bottom=-0.15,top=1.2)

###==============###==============###==============###

# ax = Axes3D(fig)
ax[:autoscale](enable = true, tight = true)

###==============###==============###==============###

fig = plt[:figure]()
fig[:set_size_inches](sx, sy, forward = true)

ax = fig[:gca](projection = "3d")
ax = fig[:add_subplot](111)

plt[:clf]()
###==============###==============###==============###

# for i in rand(1:N, 256)
for i in 1:N
    plot3D(xs = x[:, i], ys= y[:, i], zs = z[:, i], zdir = "z", lw = 0.5)
    # plot3D(x[:, i], y[:, i], color = "0.5" ,zdir="z", zs=minimum(z), lw = 0.2)
end


###==============###==============###==============###
# for i in rand(1:N, 256)
for i in 1:N
    plt[:plot](x[:, i], y[:, i], lw = 0.3)
end

###==============###==============###==============###

ax[:grid](false)

ax[:xaxis][:pane][:set_edgecolor]("white")
ax[:yaxis][:pane][:set_edgecolor]("white")
ax[:zaxis][:pane][:set_edgecolor]("white")

ax[:xaxis][:pane][:fill] = false
ax[:yaxis][:pane][:fill] = false
ax[:zaxis][:pane][:fill] = true

ax[:set_xlabel](L"\mathrm{x \times 10^4}", labelpad =-8, fontsize = fs)
ax[:set_ylabel](L"\mathrm{y \times 10^4}", labelpad =-8, fontsize = fs)
ax[:set_zlabel](L"\mathrm{z \times 10^4}", labelpad =-8, fontsize = fs)

ax[:set_xlabel](L"\mathrm{x \times 10^3}", labelpad =-8, fontsize = fs)
ax[:set_ylabel](L"\mathrm{y \times 10^3}", labelpad =-8, fontsize = fs)
ax[:set_zlabel](L"\mathrm{z \times 10^3}", labelpad =-8, fontsize = fs)

ax[:set_xlabel](L"\mathrm{x \times 10^4}", labelpad =0, fontsize = fs)
ax[:set_ylabel](L"\mathrm{y \times 10^4}", labelpad = 2, fontsize = fs)

ax[:relim]()
ax[:autoscale_view](tight = true)

fig[:subplots_adjust](top=1.1)
fig[:subplots_adjust](bottom=0)
fig[:subplots_adjust](right=0.9)
fig[:subplots_adjust](left=0)

ax[:elev] = 16
ax[:azim] = -45

###==============###==============###==============###

plt[:tick_params](which = "both", labelsize = ls, direction = "in", pad = -3)
plt[:tick_params](which = "both", labelsize = ls, direction = "in", pad = 2)

### K = 0

ax[:set_xticks](collect(linspace(-15exp10(4),0, 4)))
ax[:set_yticks](collect(linspace(-2exp10(4),8exp10(4), 6)))
ax[:set_zticks](collect(linspace(-3exp10(3),0,4)))

ax[:set_xticklabels](collect(-15:5:0))
ax[:set_yticklabels](collect(-2:2:8))
ax[:set_zticklabels](collect(-3:0))

ax[:set_zlim]([-1.2exp10(5), 0.0])
ax[:set_ylim]([-2.25exp10(5), 0.0])
###==============###==============###==============###
# 2D

plt[:tick_params](which = "both", labelsize = ls, direction = "in", pad = 1.5)

### K = 0
ax[:set_xticks]([-1exp10(4), 0 , exp10(4)])
ax[:set_yticks]([-1exp10(4), 0 , exp10(4)])

ax[:set_xticklabels](["-1", "0", "1"])
ax[:set_yticklabels](["-1", "0", "1"])

ax[:text](exp10(4), -1.2exp10(4), L"x \times 10^4", ha="center", va="center", size=fs)
ax[:text](-0.8exp10(4), exp10(4), L"y \times 10^4", ha="center", va="center", size=fs)

### K = 1.0, r = 4
ax[:set_xticks]([0 , 2exp10(4), 4exp10(4)])
ax[:set_yticks]([0 , 2exp10(4), 4exp10(4)])

ax[:set_xticklabels](["0", "2", "4"])
ax[:set_yticklabels](["0", "2", "4"])

ax[:text](3.8exp10(4), -1.1exp10(4), L"x \times 10^4", ha="center", va="center", size=fs)
ax[:text](4.2exp10(3), 3.6exp10(4), L"y \times 10^4", ha="center", va="center", size=fs)

### K = 5.0, r = 7
ax[:set_xticks]([0 , 5exp10(4), 1exp10(5)])
ax[:set_yticks]([0 , 1exp10(5), 2exp10(5)])

ax[:set_xticklabels](["0", "5", "10"])
ax[:set_yticklabels](["0", "10", "20"])

ax[:text](exp10(5), exp10(4), L"x \times 10^4", ha="center", va="center", size=fs)
ax[:text](exp10(4), 2exp10(5), L"y \times 10^4", ha="center", va="center", size=fs)

### K = 0.0, r = 2 top
ax[:set_xticks]([-4exp10(4), 0, 4exp10(4)])
ax[:set_yticks]([-4exp10(4), 0, 4exp10(4)])

ax[:set_xticklabels](["-4e+4", "0", "4e+4"])
ax[:set_yticklabels](["-4e+4", "0", "4e+4"])

ax[:text](4.2exp10(4), -4exp10(4), L"x", ha="center", va="center", size=fs)
ax[:text](-4exp10(4), 4exp10(4), L"y", ha="center", va="center", size=fs)

ax[:text](4.2exp10(4), -4exp10(4), L"x", ha="center", va="center", size=fs)
ax[:text](-4exp10(4), 4exp10(4), L"y", ha="center", va="center", size=fs)

### K = 0.025, r = 1 top
ax[:set_xticks]([0, 1exp10(5), 2exp10(5)])
ax[:set_yticks]([0, 2.5exp10(4), 5exp10(4), 7.5exp10(4)])

ax[:set_xticklabels](["0", "1e+5", "2e+5"])
ax[:set_yticklabels](["0", "2.5e+4", "5e+4", "7.5e+4"])

ax[:text](2.5exp10(5), -exp10(4), L"x", ha="center", va="center", size=fs)
ax[:text](0, 8.3exp10(4), L"y", ha="center", va="center", size=fs)

### K = 1.0, r = 5 top
ax[:set_xticks]([-5exp10(4), 0, 5exp10(4)])
ax[:set_yticks]([-5exp10(4), 0, 5exp10(4)])

ax[:set_xticklabels](["-5e+4", "0", "5e+4"])
ax[:set_yticklabels](["-5e+4", "0", "5e+4"])

ax[:text](5.3exp10(4), -6.5exp10(4), L"x", ha="center", va="center", size=fs)
ax[:text](-7.5exp10(4), 5exp10(4), L"y", ha="center", va="center", size=fs)

### K = 5.0, r = 5 top
ax[:set_xticks]([0, 5exp10(4), 1exp10(5)])
ax[:set_yticks]([0, 1exp10(5), 2exp10(5)])

ax[:set_xticklabels](["0", "5e+4", "1e+5"])
ax[:set_yticklabels](["0", "1e+5", "2e+5"])

ax[:text](1.15exp10(5), 5exp10(3), L"x", ha="center", va="center", size=fs)
ax[:text](1, 2exp10(5), L"y", ha="center", va="center", size=fs)

plt[:tight_layout]()
ax[:set_aspect]("equal")
ax[:set_aspect](1.5)
plt[:clf]()

###==============###==============###==============###

fig[:savefig]("fase_test2.eps", dpi = 100, format = "eps", bbox_inches = "tight" , pad_inches = 0.27)

fig[:savefig]("fase_test1.eps", dpi = 600, format = "eps", bbox_inches = "tight")

fig[:savefig]("fase_test1.eps", format = "eps", bbox_inches = "tight")

fig[:savefig]("2D_met_k_0.eps", format = "eps", bbox_inches = "tight")
fig[:savefig]("2D_met_k_1.eps", format = "eps", bbox_inches = "tight")
fig[:savefig]("2D_met_k_5.eps", format = "eps", bbox_inches = "tight")

fig[:savefig]("2D_top_k_0.eps", format = "eps", bbox_inches = "tight")
fig[:savefig]("2D_top_k_025.eps", format = "eps", bbox_inches = "tight")
fig[:savefig]("2D_top_k_1.eps", format = "eps", bbox_inches = "tight")
fig[:savefig]("2D_top_k_005.eps", format = "eps", bbox_inches = "tight")

fig[:savefig]("inertial_nloc_k_0.eps", format = "eps")
fig[:savefig]("inertial_nloc_k_005.eps", format = "eps")
fig[:savefig]("inertial_nloc_k_075.eps", format = "eps")

fig[:savefig]("3D_top_k_0.eps", format = "eps", bbox_inches = "tight")
fig[:savefig]("3D_top_k_0.eps", format = "eps")
fig[:savefig]("3D_top_k_0.png", format = "png", bbox_inches = "tight")

fig[:savefig]("3D_top_k_001_r9.eps", format = "eps", bbox_inches = "tight")
fig[:savefig]("3D_top_k_001_r9.eps", format = "eps")
fig[:savefig]("3D_top_k_01_r8.eps", format = "eps")
fig[:savefig]("3D_top_k_001.png", format = "png", bbox_inches = "tight")

fig[:savefig]("3D_top_k_01_r8.eps", format = "eps", bbox_inches = "tight")
fig[:savefig]("3D_top_k_01.png", format = "png", bbox_inches = "tight")

###==============###==============###==============###


###==============###==============###==============###
