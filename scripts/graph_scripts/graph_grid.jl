using Plots, CollectiveDynamics.DataAnalysis

gr()
pyplot()

### ================================== ###

N = 1024
N = 512
N = 256
N = 128
N = 64

T = 5
tau = get_times(T)

v0 = 0.1
### ================================== ###

folder_path = "$(homedir())/art_DATA/TFLOCK_DATA/DATA/data_N_$(N)"
folder_path = "$(homedir())/art_DATA/TFLOCK_DATA/L_DATA/data_N_$(N)"

# files = filter(x -> match(r"._(\d+.\d+).dat", x).captures[1] == η , readdir(folder_path))
folders = readdir(folder_path)

η_vals = unique([match(r"^\w+(\d+\.\d+)", f).captures[1] for f in folders])
T_vals = unique([match(r"._(\d+\.\d+)$", f).captures[1] for f in folders])

trays_plots = Any[]
exp_plots   = Any[]
order_plots = Any[]

for i in η_vals, j in T_vals

    data_path = string(folder_path, "/eta_", i, "_T_", j)
    println("/eta_", i, "_T_", j)

    means = Array{Float64}[]
    psi   = Array{Float64}[]

    raw_data = reinterpret(Float64, read(data_path * "/pos_1.dat"))
    pos_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))
    # pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

    raw_data = reinterpret(Float64, read(data_path * "/vel_1.dat"))
    vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))
    # vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

    x = view(pos_data, :, 1:3:3N)
    y = view(pos_data, :, 2:3:3N)
    z = view(pos_data, :, 3:3:3N)

    push!(trays_plots, plot(x, y, z, leg = false, tickfont = font(1)))

    pos_data[:, 1:3:3N] .-= mean(pos_data[:, 1:3:3N], 2)
    pos_data[:, 2:3:3N] .-= mean(pos_data[:, 2:3:3N], 2)
    pos_data[:, 3:3:3N] .-= mean(pos_data[:, 3:3:3N], 2)

    tr_pos_data = transpose(pos_data)

    means = [mean(calc_rij_3D_vect(tr_pos_data[:, i])) for i in 1:size(tr_pos_data,2)]

    psi = (1. / v0) * [norm(mean([ [vel_data[i, j], vel_data[i+1, j], vel_data[i+2, j] ] for i in 1:3:3N])) for j in 1:size(vel_data, 2)]

    push!(exp_plots, plot(tau, means, xscale = :log10, marker = :o, markersize = 1.2, leg = false, xlabel = "$(i)", ylabel = "$(j)", tickfont = font(8)))

    push!(order_plots, plot(tau, psi, xscale = :log10, marker = :o, markersize = 1.2, leg = false, xlabel = "$(i)", ylabel = "$(j)", tickfont = font(8)))

end

l = @layout [  a{0.5w} [b
    c]]

plot( trays, expansion, order, layout = l )

l = @layout grid(length(η_vals), length(T_vals))

plot(trays_plots..., layout = l,  size = (1270,820))
savefig("$(homedir())/Google\ Drive/proyecto_martin/imagenes/modelo_cvgn/trays_N_$(N).png")

plot(exp_plots..., layout = l, size = (1270,820))
savefig("$(homedir())/Google\ Drive/proyecto_martin/imagenes/modelo_cvgn/exp_N_$(N).png")

plot(order_plots..., layout = l, size = (1270,820))
savefig("$(homedir())/Google\ Drive/proyecto_martin/imagenes/modelo_cvgn/order_N_$(N).png")

gui()
