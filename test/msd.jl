using Plots, CollectiveDynamics.DataAnalysis, LaTeXStrings
using PyPlot, CollectiveDynamics.DataAnalysis, LaTeXStrings

### ================================== ###

gr()
pyplot()
gui()

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

tau = get_times(τ)

v0 = 1.0
### ================================== ###

folder = "NLOC_DATA"
folder = "NLOC_DATA_3D"
folder = "NLOC_MET_3D"
folder = "NLOC_TOP_3D"
folder = "NLOC_TOP_3D_MEAN"
folder = "TFLOCK_NLOC_DATA"
folder = "TFLOCK_DATA"

folder_path = "$(homedir())/art_DATA/$(folder)/EXP/exp_data_N_$(N)"
folder_path = "$(homedir())/exp_DATA/$(folder)/EXP/exp_data_N_$(N)"

folder_path = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)"

# files = filter(x -> match(r"._(\d+.\d+).dat", x).captures[1] == η , readdir(folder_path))
eta_folders = readdir(folder_path)
folders = readdir(folder_path * "/" * eta_folders[1])

# η_vals = [match(r"\w+\d+\w+(\d+\.\d+)", f).captures[1] for f in folders]
# all_means = Dict()
### ================================== ###

f = 6
data_path = folder_path * "/" * eta_folders[1] * "/" * folders[f]
data_path = folder_path * "/" * eta_folders[f]
data_path = folder_path * "/" * folders[f]

reps = [match(r"\w+(\d+).\w+", x).captures[1]  for x in filter(x -> ismatch(r"^pos_", x), readdir(data_path))]

r = rand(reps)
# for r in reps

means = Array{Float64}[]
psi   = Array{Float64}[]

raw_data = reinterpret(Float64, read(data_path * "/pos_$(r).dat"))
pos_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))
# pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

x = view(pos_data, :, 1:3:3N)
y = view(pos_data, :, 2:3:3N)
z = view(pos_data, :, 3:3:3N)

# trays = plot(x, y, z, leg = false, size = (600,600), aspect_ratio = :equal)
trays = plot(x, y, z, leg = false, size = (800, 800), tickfont = font(14, "serif"),
    xticks = (linspace(-2exp10(4),2exp10(4), 5), [L"-2 \times 10^4", L"-10^4", "0", L"10^4", L"2 \times 10^4" ]),
    yticks = (linspace(-2exp10(4),2exp10(4), 5), [L"-2 \times 10^4", L"-10^4", "0", L"10^4", L"2 \times 10^4" ]),
    zticks = (linspace(-2exp10(4),2exp10(4), 5), [L"-2 \times 10^4", L"-10^4", "0", L"10^4", L"2 \times 10^4" ]),
    grid = false)

savefig("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)/N_$(N)/k_0_125.png")
savefig("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)/N_$(N)/test_1.svg")

gr()
pyplot()
gui()
