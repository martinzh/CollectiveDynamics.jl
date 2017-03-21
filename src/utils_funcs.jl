## =========================== ## ## =========================== ##
## 	   Utility functions for CollectiveDynamics statistics		 ##
##	   Martin Zumaya Hernandez 						             ##
##     21 / 03 / 2017									         ##
## =========================== ## ## =========================== ##

function test()
    println("test submodule function")
end

### ============== ### ============== ### ============== ###
##           COMPUTE RELATIVE DISTANCES VECTOR            ##
### ============== ### ============== ### ============== ###
"""
    calc_rij_3D_vect(vect)
Compute distance between particles and store them in a 1-D array
# Arguments
* vect -> array of particles' positions
"""
function calc_rij_3D_vect(vect)

    n = div(length(vect),3)

    vec_rij = zeros(Float64, div(n*(n-1),2) )

    k = 1

    for i in 1:3:3n, j in (i+3):3:3n
        vec_rij[k] = norm([vect[i] - vect[j], vect[i+1] - vect[j+1], vect[i+2] - vect[j+2]])
        k += 1
    end

    return vec_rij
end

### ================================== ###

"""
    calc_rij_2D_vect(vect)
Compute distance between particles and store them in a 1-D array
# Arguments
* vect -> array of particles' positions
"""
function calc_rij_2D_vect(vect)

    n = div(length(vect), 2)

    vec_rij = zeros(Float64, div(n*(n-1),2) )

    k = 1

    for i in 1:2:2n, j in (i+2):2:2n
        vec_rij[k] = norm([vect[i] - vect[j], vect[i+1] - vect[j+1]])
        k += 1
    end

    return vec_rij
end

### ============== ### ============== ### ============== ###
##  COMPUTE POSITIONS IN CENTER OF MASS REFERENCE FRAME   ##
### ============== ### ============== ### ============== ###
"""
    calc_vect_3D_cm(vals)
Compute positions in center of mass referenece frame
"""
function calc_vect_3D_cm(vals)

    N = size(vals, 1)

    x_cm = mean(vals[1:3:N, :], 1)
    y_cm = mean(vals[2:3:N, :], 1)
    Z_cm = mean(vals[3:3:N, :], 1)

    for i in 1:size(vals,2)
        vals[1:3:N, i] -= x_cm[i]
        vals[2:3:N, i] -= y_cm[i]
        vals[3:3:N, i] -= z_cm[i]
    end
end

### ================================== ###
"""
    calc_vect_2D_cm(vals)
Compute positions in center of mass referenece frame
"""
function calc_vect_2D_cm(vals)

    N = size(vals, 1)

    x_cm = mean(vals[1:2:N, :], 1)
    y_cm = mean(vals[2:2:N, :], 1)

    for i in 1:size(vals,2)
        vals[1:2:N, i] -= x_cm[i]
        vals[2:2:N, i] -= y_cm[i]
    end
end

### ================================== ###

### ============== ### ============== ### ============== ###
##                 COMPUTE SAMPLING TIMES                 ##
### ============== ### ============== ### ============== ###
"""
    get_times(T)
Compute sampling times for 10^T time steps
"""
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

### ============== ### ============== ### ============== ###
##             FINITE DIFFERENCES DERIVATIVES             ##
### ============== ### ============== ### ============== ###
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
