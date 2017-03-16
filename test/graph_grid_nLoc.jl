using Plots

gr()
pyplot()

### ================================== ###

function calc_rij_vect(vect)

    n = div(length(vect),3)

    # println(size(vect))

    vec_rij = zeros(Float64, div(n*(n-1),2) )

    k = 1

    for i in 1:3:3n, j in (i+3):3:3n

        # vec_rij[k] = sqrt( (vect[i] - vect[j])^2 + (vect[i+1] - vect[j+1])^2 + (vect[i+2] - vect[j+2])^2)

        vec_rij[k] = norm([vect[i] - vect[j], vect[i+1] - vect[j+1], vect[i+2] - vect[j+2]])

        k += 1
    end

    # println("Finish calc rij")

    return vec_rij

end

### ================================== ###

function get_times(T)

    times = [convert(Int, exp10(i)) for i in 0:T]

    tau = Int64[]

    for i in 1:(length(times) - 1)

        if i > 1

            for t in (times[i]+1):times[i+1]

                if t % times[i] == 0 || t % times[i-1] == 0
                    push!(tau, t)
                end
            end

        else

            for t in (times[i]+1):times[i+1]

                if t % times[i] == 0
                    push!(tau, t)
                end
            end

        end

    end

    return tau
end


function deriv_interpol(vals, time)

    n = length(vals)
    der = zeros(n)

    for i in 1:3:n

        h_1 = vals[i+1] - vals[i]
        h_2 = vals[i+2] - vals[i]

        der[i] = ((2*h_1 + h_2)/(h_1*(h_1+h_2))) * vals[i] + ((h_1 + h_2)/(h_1 * h_2)) * vals[i+1] - ((h_1)/(h_2 * (h_1 + h_2))) * vals[i+2]
        der[i + 1] = -(h_2/(h_1*(h_1+h_2))) * vals[i] - ((h_1 - h_2)/(h_1 * h_2)) * vals[i+1] - (h_1/(h_2 * (h_1 + h_2))) * vals[i+2]
        der[i + 2] = (h_2/(h_1*(h_1+h_2))) * vals[i] - ((h_1 + h_2)/(h_1 * h_2)) * vals[i+1] + ((h_1 + 2*h_2)/(h_2 * (h_1 + h_2))) * vals[i+2]

    end

    return der
end

### ================================== ###

N = 1024
N = 512
N = 256
N = 128
N = 64

T = 5
T = 6
tau = get_times(T)

v0 = 0.1
### ================================== ###

folder_path = "$(homedir())/art_DATA/TFLOCK_NLOC_DATA/DATA/data_N_$(N)"

# files = filter(x -> match(r"._(\d+.\d+).dat", x).captures[1] == η , readdir(folder_path))
folders = readdir(folder_path)

η_vals = unique([match(r"^\w+(\d+\.\d+)", f).captures[1] for f in folders])

T_vals = unique([match(r"\w+_\d\.\d\w+_(\d+\.\d+)_\w+_\d$", f).captures[1] for f in readdir(folder_path * "/" * folders[1]) ] )

n_nl = 4

i = 1
# for i in 1:length(folders)

trays_plots = Any[]
exp_plots   = Any[]
order_plots = Any[]

println(folders[i])

data_path = folder_path * "/" * folders[i] * "/" * folders[i] * "_T_0.01_nl_$(n_nl)"

means = Array{Float64}[]
psi   = Array{Float64}[]

reps = unique([match(r"\w+_(\d+)\.\w+$", f).captures[1] for f in readdir(data_path)])

for r in reps

    println(r)

    raw_data = reinterpret(Float64, read(data_path * "/pos_$(r).dat"))
    pos_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))
    # pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

    raw_data = reinterpret(Float64, read(data_path * "/vel_$(r).dat"))
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

    means = [mean(calc_rij_vect(tr_pos_data[:, i])) for i in 1:size(tr_pos_data,2)]

    psi = (1. / v0) * [norm(mean([ [vel_data[i, j], vel_data[i+1, j], vel_data[i+2, j] ] for i in 1:3:3N])) for j in 1:size(vel_data, 2)]

    push!(exp_plots, plot(tau, means, xscale = :log10, marker = :o, markersize = 1.2, leg = false, xlabel = "$(r)", tickfont = font(8)))

    push!(order_plots, plot(tau, psi, xscale = :log10, marker = :o, markersize = 1.2, leg = false, xlabel = "$(r)", tickfont = font(8)))
end

plot(trays_plots...,  size = (1270,820))
savefig("$(homedir())/Google\ Drive/proyecto_martin/imagenes/modelo_cvgn/reps_trays_N_$(N)_$(folders[i])_T_0.01_nl_$(n_nl).png")

plot(exp_plots..., size = (1270,820))
savefig("$(homedir())/Google\ Drive/proyecto_martin/imagenes/modelo_cvgn/reps_exp_N_$(N)_$(folders[i])_T_0.01_nl_$(n_nl).png")

plot(order_plots..., size = (1270,820))
savefig("$(homedir())/Google\ Drive/proyecto_martin/imagenes/modelo_cvgn/reps_order_N_$(N)_$(folders[i])_T_0.01_nl_$(n_nl).png")

# end

gui()
