### ============== ### ============== ### ============== ###
### COMPUTE RELATIVE DISTANCES
### ============== ### ============== ### ============== ###

@time @parallel for i in 1:3:length(pos)
        for j in (i+3):3:length(pos)
            R_ij[div(j,3)+1,div(i,3)+1] = 0.0
            d = (pos[i]-pos[j])^2 + (pos[i+1]-pos[j+1])^2 + (pos[i+2]-pos[j+2])^2
            d <= r0^2 && d > 0.0 ? R_ij[div(j,3)+1,div(i,3)+1] = 1.0 : R_ij[div(j,3)+1,div(i,3)+1] = -1.0
        end
end

### ============== ### ============== ### ============== ###
### COMPUTE SHORT AND LONG RANGE INTERACTIONS
### ============== ### ============== ### ============== ###

function compute_metric_interactions(vel::SharedArray,v_r::SharedArray,v_n::SharedArray)

    for id in first(localindexes(vel)):3:last(localindexes(vel))

        i = div(id, 3)
        # print(i+1,"|\t")

        v_r[3i+1] = 0.0
        v_r[3i+2] = 0.0
        v_r[3i+3] = 0.0

        v_n[3i+1] = 0.0
        v_n[3i+2] = 0.0
        v_n[3i+3] = 0.0

        sh_n = find(x -> x > zero(Float64), Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])
        k_sh = length(sh_n)
        # k_sh = sum(sh_n)

        ln_n = find(x -> x < zero(Float64), Symmetric(R_ij, :L)[(i*N)+1:(i+1)*N])

        # short-range
        if isempty(sh_n) == false
            for j in sh_n
                # print(j,"\t")
                v_r[3i+1] += vel[3(j-1)+1] / k_sh
                v_r[3i+2] += vel[3(j-1)+2] / k_sh
                v_r[3i+3] += vel[3(j-1)+3] / k_sh
            end
        end

        # println()
        k_ln = rand(κ_dist)
        # println(k_ln)

        # possible long range
        # if isempty(ln_n) == false && k_ln != 0.0
        if k_ln != 0.0
            for j in rand(ln_n, 3)
                # print(j,"\t")
                v_n[3i+1] += vel[3(j-1)+1] / k_ln
                v_n[3i+2] += vel[3(j-1)+2] / k_ln
                v_n[3i+3] += vel[3(j-1)+3] / k_ln
            end
        end

        # println()

    end
end

### ============== ### ============== ### ============== ###
### UPDATE PARTICLE'S POSITIONS AND VELOCITIES
### ============== ### ============== ### ============== ###

function update_particles(pos::SharedArray, vel::SharedArray,v_r::SharedArray,v_n::SharedArray)

    for id in first(localindexes(vel)):3:last(localindexes(vel))

        i = div(id, 3)
        # print(i+1,"|\t")

        signal = ω * [v_r[3i+1] , v_r[3i+2], v_r[3i+3]] + (1.0 - ω) * [v_n[3i+1] , v_n[3i+2], v_n[3i+3]]
        signal_angle = 0.0

        p_vel = [vel[3i+1] , vel[3i+2], vel[3i+3]]

        norm(signal) != zero(Float64) ? signal_angle = acos(dot(p_vel, signal) / norm(signal)) : signal_angle = 0.0

        q_r = Quaternion(zeros(Float64, 3))

        if norm(signal) != zero(Float64)
            q_r = qrotation(cross(p_vel, signal), signal_angle + η * (2.0 * rand() * pi - pi)) * Quaternion(p_vel)
        else
            noise = randn(3)
            q_r = qrotation(cross(p_vel, noise), η * acos(dot(normalize(noise), p_vel)) ) * Quaternion(p_vel)
        end

        u_vel = normalize([q_r.v1, q_r.v2, q_r.v3])

        vel[3i+1] = u_vel[1]
        vel[3i+2] = u_vel[2]
        vel[3i+3] = u_vel[3]

        pos[3i+1] += u_vel[1]
        pos[3i+2] += u_vel[2]
        pos[3i+3] += u_vel[3]

    end
end
