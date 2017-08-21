using PyPlot, CollectiveDynamics.DataAnalysis

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
folder = "NLOC_MET_2D"

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

f = 9
data_path = folder_path * "/" * eta_folders[f]

# data_path = folder_path * "/" * eta_folders[1] * "/" * folders[f]
# data_path = folder_path * "/" * folders[f]

reps = [match(r"\w+(\d+).\w+", x).captures[1]  for x in filter(x -> ismatch(r"^pos_", x), readdir(data_path))]

r = rand(reps)
# for r in reps

means = Array{Float64}[]
psi   = Array{Float64}[]

raw_data = reinterpret(Float64, read(data_path * "/pos_$(r).dat"))

pos_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))
pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

file = open("pos_3D_k0.csv", "w+")

for i in 1:size(pos_data, 2)
    println(i)
    println(file, repr(pos_data[:,i])[2:end-1])
end

close(file)

pos_data = transpose(reshape(raw_data, 2N, div(length(raw_data), 2N)))
# pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

### ================================== ###

x = view(pos_data, :, 1:3:3N)
y = view(pos_data, :, 2:3:3N)
z = view(pos_data, :, 3:3:3N)

# trays = plot(x, y, z, leg = false, size = (600,600), aspect_ratio = :equal)
# trays = plot(x, y, z, leg = false, size = (800,800), tickfont = font(1), aspect_ratio = :equal)

for i in rand(1:N, 256)
    plot3D(xs = x[:,i], ys = y[:,i], zs = z[:,i], zdir = "z", lw = 0.8)
end

plt[:ticklabel_format](style = "sci", scilimits = (2,4), useOffset = false)

plt[:NullFormatter]()

plt[:clf]()

rand(1:N, 100)

savefig("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)/N_$(N)/k_0_5.png")

### ================================== ###

x = view(pos_data, :, 1:2:2N)
y = view(pos_data, :, 2:2:2N)

# trays = plot(x, y, z, leg = false, size = (600,600), aspect_ratio = :equal)
# trays = plot(x, y, z, leg = false, size = (800,800), tickfont = font(1), aspect_ratio = :equal)

for i in rand(1:N, 256)
    plt[:plot](x[:,i], y[:,i], lw = 0.8)
end

plt[:ticklabel_format](style = "sci", scilimits = (2,4), useOffset = false)

plt[:NullFormatter]()

plt[:clf]()

rand(1:N, 100)

savefig("/Users/mzumaya/Google Drive/proyecto_martin/graphs_p_mod/$(folder)/N_$(N)/k_0_5.png")
