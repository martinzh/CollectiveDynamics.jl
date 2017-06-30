ENV["PLOTS_USE_ATOM_PLOTPLANE"] = "false"

using PyPlot, LaTeXStrings, CollectiveDynamics.DataAnalysis
using PyPlot

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
folder = "NLOC_MET_3D"
folder = "NLOC_TOP_3D"
folder = "NLOC_TOP_3D_MEAN"
folder = "TFLOCK_NLOC_DATA"
folder = "TFLOCK_DATA"

folder_path = "$(homedir())/art_DATA/$(folder)/EXP/exp_data_N_$(N)"
folder_path = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)"

eta_folders = readdir(folder_path)

f = 10
data_path = folder_path * "/" * eta_folders[f]

reps = [match(r"\w+(\d+).\w+", x).captures[1]  for x in filter(x -> ismatch(r"^pos_", x), readdir(data_path))]

r = rand(reps)

raw_data = reinterpret(Float64, read(data_path * "/pos_$(r).dat"))
pos_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))

x = view(pos_data, :, 1:3:3N)
y = view(pos_data, :, 2:3:3N)
z = view(pos_data, :, 3:3:3N)

###==============###==============###==============###

plt[:rc]("font", family="serif")
plt[:rc]("text", usetex=true)

# plt[:rc]("figure.autolayout")
# plt[:rcParams]["figure.autolayout"] = true

###==============###==============###==============###
fs = 14
ls = 10

sx = sy = 3

###==============###==============###==============###

plt[:ticklabel_format](style = "sci", scilimits = (2,4))
plt[:ticklabel_format](style = "plain")

fig = plt[:figure](num = 1, dpi = 100, facecolor="w", edgecolor="k")
fig[:set_size_inches](sx, sy, forward = true)

# fig[:subplots_adjust](bottom=-0.15,top=1.2)

###==============###==============###==============###

# ax = Axes3D(fig)
ax[:autoscale](enable = true, tight = true)

###==============###==============###==============###
fig = plt[:figure]()
fig[:set_size_inches](sx, sy, forward = true)

ax = fig[:gca](projection = "3d")

###==============###==============###==============###

for i in rand(1:N, 256)
    plot3D(xs = x[:, i], ys= y[:, i], zs = z[:, i], zdir = "z", lw = 0.5)
    # plot3D(x[:, i], y[:, i], color = "0.5" ,zdir="z", zs=minimum(z), lw = 0.2)
end

###==============###==============###==============###

ax[:grid](false)

ax[:xaxis][:pane][:set_edgecolor]("white")
ax[:yaxis][:pane][:set_edgecolor]("white")
ax[:zaxis][:pane][:set_edgecolor]("white")

ax[:xaxis][:pane][:fill] = false
ax[:yaxis][:pane][:fill] = false
ax[:zaxis][:pane][:fill] = true

ax[:set_xlabel](L"\mathrm{x}", labelpad =0)
ax[:set_ylabel](L"\mathrm{y}", labelpad =0)
ax[:set_zlabel](L"\mathrm{z}", labelpad =0)

ax[:relim]()
ax[:autoscale_view](tight = true)
fig[:subplots_adjust](top=1.05)
###==============###==============###==============###

plt[:tick_params](which = "both", labelsize = ls, direction = "in", pad = 1.5)

### K = 0

ax[:set_xticks](collect(linspace(-2exp10(4),2exp10(4),5)))
ax[:set_yticks](collect(linspace(0, exp10(5),4)))
ax[:set_zticks](collect(linspace(-1exp10(4),2exp10(4),5)))

ax[:set_xticklabels](collect(-2:1:2))
ax[:set_yticklabels](collect(-2:1:2))
ax[:set_zticklabels](collect(-2:1:2))

ax[:set_xticklabels](["-1e+4", "0", "1e+4"])
ax[:set_yticklabels](["-1e+4", "0", "1e+4"])
ax[:set_zticklabels](["-1e+4", "0", "1e+4"])

ax[:set_xticks]([0, 5exp10(4), 1exp10(5)])
ax[:set_yticks]([0, -5exp10(4), -1.5exp10(5)])
ax[:set_zticks]([0, -4exp10(4), -8exp10(4), -1.2exp10(5)])

ax[:set_xticklabels](["-1.5e+5", "-1e+5", "-5e+4", "0"])
ax[:set_yticklabels](["0", "4e+4", "8e+4", "1.2e+5"])
ax[:set_zticklabels](["0", "5e+4", "1e+5"])

ax[:set_xticklabels](["-5e+4", "0", "5e+4"])
ax[:set_yticklabels](["0", "-5e+4", "-1.5e+5"])
ax[:set_zticklabels](["-6e+4", "-4e+4", "-2e+4", "0"])

ax[:set_zlim]([-1.2exp10(5), 0.0])
ax[:set_ylim]([-2.25exp10(5), 0.0])
###==============###==============###==============###

fig[:savefig]("fase_test2.eps", dpi = 100, format = "eps", bbox_inches = "tight" , pad_inches = 0.27)

fig[:savefig]("fase_test1.eps", dpi = 600, format = "eps", bbox_inches = "tight")

fig[:savefig]("fase_test1.eps", format = "eps", bbox_inches = "tight")

###==============###==============###==============###

plt[:clf]()

###==============###==============###==============###
