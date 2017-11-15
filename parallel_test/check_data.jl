using Plots

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

folder = "NLOC_P_MET_3D"
folder = "NLOC_P_TOP_3D"

folder_path = "$(homedir())/art_DATA/$(folder)/DATA/data_N_$(N)"
folder_path = "$(homedir())/art_DATA/$(folder)/EXP/exp_data_N_$(N)"

eta_folders = readdir(folder_path)

folders = readdir(folder_path * "/" * eta_folders[1])

### ================================== ###

f = 1
data_path = folder_path * "/" * eta_folders[f]

reps = [match(r"\w+(\d+).\w+", x).captures[1]  for x in filter(x -> ismatch(r"^pos_", x), readdir(data_path))]

r = rand(reps)

raw_data = reinterpret(Float64, read(data_path * "/pos_$(r).dat"))

# 3D
pos_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))

x = view(pos_data, :, 1:3:3N)
y = view(pos_data, :, 2:3:3N)
z = view(pos_data, :, 3:3:3N)

plot(x, y, z, leg = false, size = (800,800), tickfont = font(1), aspect_ratio = :equal)
