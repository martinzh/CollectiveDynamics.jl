using PyPlot
ENV["PLOTS_USE_ATOM_PLOTPLANE"] = "false"

###==============###==============###==============###

met_vals = readcsv("/Users/mzumaya/Google Drive/proyecto_martin/3D_met_order.csv")
top_vals = readcsv("/Users/mzumaya/Google Drive/proyecto_martin/3D_top_order.csv")

###==============###==============###==============###

plt[:rc]("font", family="serif")

fig = plt[:figure](num = 1, dpi = 300, facecolor="w", edgecolor="k")
fig[:set_size_inches](2.27, 1.4, forward = true)

ax = fig[:add_subplot](111)


###==============###==============###==============###

ax[:plot](met_vals[:, 1], met_vals[:, 2], "-o", top_vals[:, 1], top_vals[:, 2], "-^", ms = 4, lw = 0.5)

#7C98AB, #F6883D
plt[:xscale]("log")
plt[:yscale]("linear")

plt[:xlim](0.0, 1.55)

plt[:tick_params](labelsize = 10, direction = "in")
plt[:tick_params](which = "minor", direction = "in")
plt[:tight_layout]()

###==============###==============###==============###

# plt[:savefig]("fase_t.eps", dpi = "figure", format = "eps", bbox_inches = "tight" , pad_inches = 0.1)

fig[:savefig]("fase_test1.eps", dpi = 300, format = "eps", bbox_inches = "tight" , pad_inches = 0.1)

###==============###==============###==============###

plt[:clf]()

###==============###==============###==============###
