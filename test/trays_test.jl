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
folder = "NLOC_MET_3D"
folder = "NLOC_TOP_3D"
folder = "NLOC_TOP_3D_MEAN"
folder = "TFLOCK_NLOC_DATA"
folder = "TFLOCK_DATA"

folder_path = "$(homedir())/art_DATA/$(folder)/EXP/exp_data_N_$(N)"
folder_path = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)"

eta_folders = readdir(folder_path)

f = 12
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

plt[:rc]("figure.autolayout")

plt[:rcParams]["figure.autolayout"] = true

###==============###==============###==============###
fs = 14
ls = 8

sx = 2.5
sy = 2.5

###==============###==============###==============###

plt[:ticklabel_format](style = "sci", scilimits = (2,4))
plt[:ticklabel_format](style = "plain")

fig = plt[:figure](num = 1, dpi = 300, facecolor="w", edgecolor="k", autolayout = true)
fig[:set_size_inches](sx, sy, forward = true)

fig[:subplots_adjust](bottom=-0.15,top=1.2)

###==============###==============###==============###

ax = Axes3D(fig)
ax[:autoscale](enable = true)

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

plt[:zoom](0.8)

ax[:gcf]()[:subplots_adjust](bottom=0.35)
###==============###==============###==============###

plt[:tick_params](which = "both", labelsize = ls, direction = "in", pad = -3)

### K = 0

ax[:set_xticks](collect(linspace(-2exp10(4),2exp10(4),5)))
ax[:set_xticklabels](collect(-2:1:2))

ax[:set_yticks](collect(linspace(-2exp10(4),2exp10(4),5)))
ax[:set_yticklabels](collect(-2:1:2))

ax[:set_zticks](collect(linspace(-1exp10(4),2exp10(4),5)))
ax[:set_zticklabels](collect(-2:1:2))

### K = 0

ax[:set_yticks](collect(linspace(-2exp10(4),2exp10(4),5)))
ax[:set_yticklabels](collect(-2:1:2))

ax[:set_zticks](collect(linspace(5exp10(4),25exp10(4),5)))
ax[:set_zticklabels](["5e5","10e5","15e5","25e5"])

# ax[:set_xlabel](L"\mathrm{x}\;(10^4)")
# ax[:set_ylabel](L"\mathrm{y}\;(10^4)")
# ax[:set_zlabel](L"\mathrm{z}\;(10^4)")

###==============###==============###==============###

fig[:savefig]("fase_test1.eps", dpi = 300, format = "eps", bbox_inches = "tight" , pad_inches = 0.01)

###==============###==============###==============###

plt[:clf]()

###==============###==============###==============###
