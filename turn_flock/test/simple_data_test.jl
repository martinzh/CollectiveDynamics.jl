using Plots; gr()

### ================================== ###

function calc_rij_vect(vect)

    n = div(length(vect),3)

    # println(size(vect))

    vec_rij = zeros(Float64, div(n*(n-1),2) )

    k = 1

    for i in 1:3:3n, j in (i+3):3:3n

        vec_rij[k] = sqrt( (vect[i] - vect[j])^2 + (vect[i+1] - vect[j+1])^2 + (vect[i+2] - vect[j+2])^2)

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
N = 128
N = 64

T = 4

η = "0.1"
η = "15.0"
η = "60.0"

tau = get_times(T)
### ================================== ###

folder_path = "$(homedir())/Code/art_DATA/cvgn_data/data_N_$(N)/data_N_$(N)_eta_$(η)"
# folder_path = "$(homedir())/Code/art_DATA/cvgn_data/data_N_$(N)/data_N_$(N)_eta_$(η)"

folder_path = "$(homedir())/GitRepos/CollectiveDynamics.jl/turn_flock/TFLOCK_DATA/data_N_$(N)/data_N_$(N)_eta_$(η)"

# files = filter(x -> match(r"._(\d+.\d+).dat", x).captures[1] == η , readdir(folder_path))
files = readdir(folder_path)

### ================================== ###

# psi   = Array{Float64}[]

all_means = Dict()
means = Array{Float64}[]
psi   = Array{Float64}[]

r = 1
raw_data = reinterpret(Float64, read(folder_path * "/pos_$(r).dat"))
pos_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))
# pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

raw_data = reinterpret(Float64, read(folder_path * "/vel_$(r).dat"))
vel_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))
# vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

raw_data = reinterpret(Float64, read(folder_path * "/spin_$(r).dat"))
spin_data = transpose(reshape(raw_data, 3N, div(length(raw_data), 3N)))
# spin_data = reshape(raw_data, 3N, div(length(raw_data), 3N))

find(x -> isnan(x), pos_data[:, 1])
find(x -> isnan(x), vel_data[:, 1])
find(x -> isnan(x), spin_data[:, 1])

x = view(pos_data, :, 1:3:3N)
y = view(pos_data, :, 2:3:3N)
z = view(pos_data, :, 3:3:3N)

gui()
plot(x, y, z, leg = false)

x_mean = mean(pos_data[:, 1:3:3N], 2)

pos_data[:, 1:3:3N] = pos_data[:, 1:3:3N] .- mean(pos_data[:, 1:3:3N], 2)
pos_data[:, 2:3:3N] = pos_data[:, 2:3:3N] .- mean(pos_data[:, 2:3:3N], 2)
pos_data[:, 3:3:3N] = pos_data[:, 3:3:3N] .- mean(pos_data[:, 3:3:3N], 2)

tr_pos_data = transpose(pos_data)

means = [mean(calc_rij_vect(tr_pos_data[:, i])) for i in 1:size(tr_pos_data,2)]

tr_vel_data = transpose(vel_data)

psi = [norm(mean([ [tr_vel_data[i, j], tr_vel_data[i+1, j], tr_vel_data[i+2, j] ] for i in 1:3:3N])) for j in 1:size(tr_vel_data, 2)]

psi = [norm(mean([[tr_vel_data[i, j], tr_vel_data[i+1, j], tr_vel_data[i+2, j]] for i in 1:3:3N])) for j in 1:size(tr_vel_data, 2)]

size(psi)
println(j)

plot(collect(1:convert(Int, exp10(4))), means, xscale = :log10, yscale = :log10, marker = :o, markersize = 0.5)
plot(collect(1:convert(Int, exp10(4))), psi, xscale = :log10, marker = :o)
