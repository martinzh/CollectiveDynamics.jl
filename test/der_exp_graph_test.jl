
using Plots, CollectiveDynamics.DataAnalysis, Polynomials
# using Polynomials

### ================================== ###
function deriv_simple(Y, X)

    N = length(Y)

    δY = zeros(N-1)

    for i in 1:(N-1)

        h = X[i+1] - X[i]
        δY[i] = (Y[i+1] - Y[i]) / h

    end

    return δY
end
### ================================== ###
function deriv_sym(Y, X)

    N = length(Y)

    δY = zeros(N)

    # for i in 2:(N-1)
    for i in 1:N

        # h = X[i+1] - X[i-1]
        # δY[i] = (Y[i+1] - Y[i-1]) / (2 * h)

        if i == 1
            h = X[i+1] - X[i]
            δY[i] = (Y[i+1] - Y[i]) / (h)
        elseif i == N
            h = X[i] - X[i-1]
            δY[i] = (Y[i] - Y[i-1]) / (h)
        else
            h = X[i+1] - X[i-1]
            δY[i] = (Y[i+1] - Y[i-1]) / (2 * h)
        end
    end

    return δY
end

### ================================== ###
function deriv_interpol(fX, X)

    n = length(fX)

    δfX = zeros(n)

    for i in 1:3:(n-3)

        h_1 = X[i+1] - X[i]
        h_2 = X[i+2] - X[i+1]

        δfX[i] = ((2*h_1 + h_2)/(h_1*(h_1+h_2))) * fX[i] + ((h_1 + h_2)/(h_1 * h_2)) * fX[i+1] - ((h_1)/(h_2 * (h_1 + h_2))) * fX[i+2]
        δfX[i + 1] = -(h_2/(h_1*(h_1+h_2))) * fX[i] - ((h_1 - h_2)/(h_1 * h_2)) * fX[i+1] - (h_1/(h_2 * (h_1 + h_2))) * fX[i+2]
        δfX[i + 2] = (h_2/(h_1*(h_1+h_2))) * fX[i] - ((h_1 + h_2)/(h_1 * h_2)) * fX[i+1] + ((h_1 + 2*h_2)/(h_2 * (h_1 + h_2))) * fX[i+2]

    end
    return δfX
end
### ================================== ###
function polyfit_test(x, y, n)
  A = [ float(x[i])^p for i = 1:length(x), p = 0:n ]
  A \ y
end
### ================================== ###


gr()
pyplot()

N = 1024
N = 256

τ = 7
τ = 6

times = get_times(τ)

folder = "NLOC_DATA"
folder = "NLOC_TOP_3D"
folder = "NLOC_TOP_3D_MEAN"

folder_path = "$(homedir())/art_DATA/$(folder)/EXP/exp_data_N_$(N)"

eta_folders = readdir(folder_path)

order_files = filter( x -> ismatch(r"^order.", x), readdir(folder_path))
exp_files = filter( x -> ismatch(r"^exp.", x), readdir(folder_path))

# k_vals = [ match(r"(\d+\.\d+)\w+\d+\.\d+.dat$", x).captures[1] for x in order_files]
k_vals = [ match(r"(\d+\.\d+).dat$", x).captures[1] for x in order_files]
k_vals = [ match(r"^\w+(\d+\.\d+)_.", x).captures[1] for x in order_files]

means = zeros(length(times), length(order_files))
orders = zeros(length(times), length(order_files))
std_means = zeros(length(times), length(order_files))

# i = 1
for i in 1:length(order_files)

    raw_data = reinterpret(Float64, read(folder_path * "/" * exp_files[i]))
    exp_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    raw_data = reinterpret(Float64, read(folder_path * "/" * order_files[i]))
    order_data = reshape(raw_data, length(times), div(length(raw_data), length(times)))

    means[:, i] = mean(exp_data, 2)
    std_means[:, i] = std(exp_data, 2)
    orders[:, i] = mean(order_data, 2)

end

order_p = plot(times, orders, lab = reshape(k_vals, 1, length(k_vals)), xscale = :log10)
exp_p   = plot(times, means, lab = reshape(k_vals, 1, length(k_vals)), xscale = :log10, yscale = :log10)
exp_p   = plot(times, means, lab = reshape(k_vals, 1, length(k_vals)))
exp_p   = plot(times, means, yerror = std_means, leg   = false, xscale = :log10, yscale = :log10, size = (1024,720))

plot(exp_p, order_p, link = :x, layout = @layout [a b])

gui()

plot(times, means[:,end], xscale = :log10, yscale = :log10)

Δtau = diff(times)
Δs = unique(Δtau)

all_derivs = zeros(length(Δtau), size(means, 2))
all_derivs = zeros(542, size(means, 2))
all_derivs = zeros(453, size(means, 2))

for k in 1:size(means, 2)

    deriv = Array{Float64}[]

    for j in 1:length(Δs)
        filt_x = [times[i] for i in findin(Δtau, Δs[j])]
        filt_y = [means[i, k] for i in findin(Δtau, Δs[j])]

        # deriv = vcat(deriv, deriv_simple(filt_y, filt_x) )
        # deriv = vcat(deriv, deriv_simple(broadcast(x -> log10(x), filt_y), broadcast(x -> log10(x), filt_x)) )
        deriv = vcat(deriv, deriv_sym(broadcast(x -> log10(x), filt_y), broadcast(x -> log10(x), filt_x)) )
    end

    all_derivs[:, k] = deriv

    end

derivs_p = plot(times[1:length(all_derivs[:,1])], all_derivs, xscale = :log10, lab = reshape(k_vals, 1, length(k_vals)))

plot(exp_p, derivs_p, link = :x, layout = @layout [a b])

gui()

derivs_p = scatter(times[1:length(all_derivs[:,1])], all_derivs[:,end], xscale = :log10, alpha = 0.3, leg = false)
exp_p = plot(times, means[:,end], lab = reshape(k_vals, 1, length(k_vals)), xscale = :log10, yscale = :log10)

plot(exp_p, derivs_p, link = :x, layout = @layout [a b])

d_fit = polyfit(broadcast(x -> log10(x), times[1:length(all_derivs[:,1])]), all_derivs[:, end])

d_fit = polyfit_test(broadcast(x -> log10(x), times[1:length(all_derivs[:,1])]), all_derivs[:, end], 12)

d_fit = polyfit_test(times[1:length(all_derivs[:,1])], all_derivs[:, end], 12)

fit_vals = [polyval(Poly(d_fit), x) for x in broadcast(x -> log10(x), times[1:length(all_derivs[:,1])])]

fit_vals = [polyval(Poly(d_fit), x) for x in times[1:length(all_derivs[:,1])]]

plot(times[1:length(all_derivs[:,1])], fit_vals, xscale = :log10 )

derivs_p = scatter(times[1:length(all_derivs[:,1])], all_derivs[:,end], xscale = :log10, alpha = 0.3, leg = false)
plot!(derivs_p, times[1:length(all_derivs[:,1])], fit_vals, xscale = :log10 )

### ================================== ###
all_derivs = zeros(length(times), size(means, 2))

x_vals = broadcast(x -> log10(x), times)
x_vals = times

j = 5
for j in 1:size(means, 2)
    y_vals = broadcast(x -> log10(x), means[:, j])
    # y_vals = means[:, j]

    # ajusta valores
    p_fit = polyfit(x_vals, y_vals, 20)

    # calcula derivada de ajuste
    all_derivs[:, j] = polyval(polyder(p_fit), x_vals)
end

m_p = scatter(x_vals, means, leg = false, alpha = 0.3)
m_p = plot(times, means, leg = false, xscale = :log10, yscale = :log10)
plot!(m_p, x_vals, all_derivs)

e_lag = 50
# d_p = plot(times, all_derivs[2:end, :], leg = false, xscale = :log10)
d_p = plot(x_vals[2:end-e_lag], all_derivs[2:end-e_lag, :], leg = false)

[times[findmax(all_derivs[10:end-e_lag,i])[2]] for i in 1:size(all_derivs, 2)]



d2_p = plot(x_vals, polyval(polyder(polyder(p_fit)), x_vals))

[exp10(y) for y in sort([real(x) for x in roots(polyder(polyder(p_fit)))])]

plot(m_p, d_p, layout = @layout [a b])

# plot(m_p, d_p, d2_p, layout = @layout [a b c])
# plot(d_p, d2_p, layout = @layout [a; b])

gui()
