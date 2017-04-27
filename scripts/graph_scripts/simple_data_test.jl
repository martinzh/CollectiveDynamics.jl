using Plots, CollectiveDynamics.DataAnalysis

### ================================== ###

gr()
pyplot()
gui()

### ================================== ###

N = 1024
N = 512
N = 256
N = 128
N = 64

τ = 5
τ = 4
τ = 3

tau = get_times(τ)

v0 = 1.0
### ================================== ###

folder_path = "$(homedir())/art_DATA/NLOC_DATA_3D/DATA/data_N_$(N)"
folder_path = "$(homedir())/art_DATA/NLOC_TOP_3D/DATA/data_N_$(N)"

# files = filter(x -> match(r"._(\d+.\d+).dat", x).captures[1] == η , readdir(folder_path))
folders = readdir(folder_path)

# η_vals = [match(r"\w+\d+\w+(\d+\.\d+)", f).captures[1] for f in folders]
# all_means = Dict()
### ================================== ###

f = 1
data_path = folder_path * "/" * folders[f]

reps = [match(r"\w+(\d+).\w+", x).captures[1]  for x in filter(x -> ismatch(r"^pos_", x), readdir(data_path))]

r = 1
# for r in reps

means = Array{Float64}[]
psi   = Array{Float64}[]

raw_data = reinterpret(Float64, read(data_path * "/pos_$(r).dat"))
pos_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))
# pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

x = view(pos_data, :, 1:3:3N)
y = view(pos_data, :, 2:3:3N)
z = view(pos_data, :, 3:3:3N)

trays = plot(x, y, z, leg = false)
gui()

raw_data = reinterpret(Float64, read(data_path * "/vel_$(r).dat"))
vel_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))
# vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

pos_data[:, 1:3:3N] = pos_data[:, 1:3:3N] .- mean(pos_data[:, 1:3:3N], 2)
pos_data[:, 2:3:3N] = pos_data[:, 2:3:3N] .- mean(pos_data[:, 2:3:3N], 2)
pos_data[:, 3:3:3N] = pos_data[:, 3:3:3N] .- mean(pos_data[:, 3:3:3N], 2)

tr_pos_data = transpose(pos_data)

means = [mean(calc_rij_3D_vect(tr_pos_data[:, i])) for i in 1:size(tr_pos_data,2)]

tr_vel_data = transpose(vel_data)

psi = (1. / v0) * [norm(mean([ [tr_vel_data[i, j], tr_vel_data[i+1, j], tr_vel_data[i+2, j] ] for i in 1:3:3N])) for j in 1:size(tr_vel_data, 2)]

expansion = plot(tau, means, xscale = :log10, yscale = :log10, marker = :o, markersize = 2.0, leg = false, title = "expansion", titlefont = font(10))

order = plot(tau, psi, xscale = :log10, marker = :o, markersize = 2.0, leg = false, title = "orden", titlefont = font(10))

l = @layout [  a{0.5w} [b
    c]]

plot( trays, expansion, order, layout = l )
