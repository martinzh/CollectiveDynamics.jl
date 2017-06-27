ENV["PLOTS_USE_ATOM_PLOTPLANE"] = "false"

using PyPlot
using LaTeXStrings

###==============###==============###==============###

# met_vals = readcsv("/Users/mzumaya/Google Drive/proyecto_martin/3D_met_order.csv")
# top_vals = readcsv("/Users/mzumaya/Google Drive/proyecto_martin/3D_top_order.csv")

met_vals = readcsv("/Users/mzumaya/GitRepos/CollectiveDynamics.jl/3D_met_order.csv")
top_vals = readcsv("/Users/mzumaya/GitRepos/CollectiveDynamics.jl/3D_top_order.csv")

###==============###==============###==============###

plt[:rc]("font", family="serif")
plt[:rc]("text", usetex=true)

fs = 14
ls = 14

###==============###==============###==============###

fig = plt[:figure](num = 1, dpi = 300, facecolor="w", edgecolor="k")
fig[:set_size_inches](2.4, 1.92, forward = true)
fig[:set_size_inches](2, 2, forward = true)

ax = fig[:add_subplot](111)

###==============###==============###==============###

ax[:plot](met_vals[:, 1], met_vals[:, 2], "-o", color = "#7C98AB", ms = 2, lw = 0.5)
ax[:plot](top_vals[:, 1], top_vals[:, 2], "-^", color = "#F6883D", ms = 2, lw = 0.5)

#7C98AB, #F6883D
plt[:xscale]("log")
# plt[:xscale]("linear")

plt[:xlim](0.0, 1.55)
plt[:yticks](collect(linspace(0,1,3)))

ax[:text](1.1, 0.1, L"\kappa", ha="center", va="center", size=fs)
ax[:text](0.0025, 0.9, L"\Psi(\kappa)", ha="center", va="center", size=fs)

plt[:tick_params](which = "both", labelsize = ls, direction = "in", pad = 4)
plt[:tight_layout]()

###==============###==============###==============###

# plt[:savefig]("fase_t.eps", dpi = "figure", format = "eps", bbox_inches = "tight" , pad_inches = 0.1)

fig[:savefig]("fase_test1.eps", dpi = 300, format = "eps", bbox_inches = "tight" , pad_inches = 0.1)

###==============###==============###==============###

plt[:clf]()

###==============###==============###==============###
