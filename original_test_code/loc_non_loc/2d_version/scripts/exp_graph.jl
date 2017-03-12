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
### ================================== ###

using Plots

gr()
pyplot()

N = 1024

τ = 7

times = get_times(τ)

folder_path = "$(homedir())/art_DATA/NLOC_DATA/EXP/exp_data_N_$(N)"

order_files = filter( x -> ismatch(r"order_k_\d+\.\d+.dat", x), readdir(folder_path))
exp_files = filter( x -> ismatch(r"exp_k_\d+\.\d+.dat", x), readdir(folder_path))

# i = 2

means = zeros(length(times), length(order_files))
orders = zeros(length(times), length(order_files))

i = 1
for i in 1:length(order_files)

    raw_data = reinterpret(Float64, read(folder_path * "/" * exp_files[i]))
    exp_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    raw_data = reinterpret(Float64, read(folder_path * "/" * order_files[i]))
    order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    means[:, i] = mean(exp_data, 2)
    orders[:, i] = mean(order_data, 2)

end


std(exp_data, 2)

order_p = plot(times, orders, leg = false, xscale = :log10)
exp_p   = plot(times, means, leg   = false, xscale = :log10, yscale = :log10)
exp_p   = plot(times, means, yerror = std(exp_data,2), leg   = false, xscale = :log10, yscale = :log10, size = (1024,720))

plot(exp_p, order_p, link = :x, layout = @layout [a ;b])

gui()
