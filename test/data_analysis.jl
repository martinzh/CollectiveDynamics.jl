using StatPlots

N = 1024
κ = 5

folder_path = "DATA/data_N_$(N)_kappa_$(κ)"

r = 5

raw_data = reinterpret(Float64,read(folder_path * "/pos_$(r).dat"))

pos_data = reshape(raw_data, 2N, div(length(raw_data), 2N))

raw_data = reinterpret(Float64,read(folder_path * "/vel_$(r).dat"))

vel_data = reshape(raw_data, 2N, div(length(raw_data), 2N))

raw_data = reinterpret(Float64,read(folder_path * "/net_$(r).dat"))

int_net = reshape(raw_data, N, N)

anim = Animation()

for t in 1:size(pos_data, 2)

    println(t)

    x = [pos_data[i, t] for i in 1:2:2N]
    y = [pos_data[i+1, t] for i in 1:2:2N]

    pts = vec(P2[(x[i],y[i]) for i in 1:N])

    # overlapping graphs
    # quiver!(pts, quiver = ([vel[i][1] for i in 1:N], [vel[i][2] for i in 1:N]))

    # overwriting graphs
    quiver(pts, quiver = ([vel_data[i, t] for i in 1:2:2N], [vel_data[i+1, t] for i in 1:2:2N]))

    frame(anim)

end


gif(anim, "test.gif", fps = 24)
