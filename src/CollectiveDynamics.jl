module CollectiveDynamics

export Flock, make_dir, make_IN, init_pos_vel, evol

###========================================###
###  FLOCK TYPE DEFINITION TO MANAGE ARGS  ###
###========================================###

type Flock
    n::Int64
    pos::Array{Float64,1}
    vels::Array{Float64,1}
    VN::Array{Float64,1}
    VR::Array{Float64,1}
    Nij::Array{Int64,1}
    poski::Array{Int64,1}

    Flock(n::Int64,m::Int64) = new(n,zeros(Float64,2n),zeros(Float64,2n),zeros(Float64,2n),zeros(Float64,2n),zeros(Int64,n+m),zeros(Int64,n))
end

###========================================###

function make_dir(path::ASCIIString)
    try
        run(`mkdir $path`)
    catch y
        println(typeof(y))
    end
end


###========================================###
###   RANDOM NUMBERS GENS (NOISE & INIT)   ###
###========================================###

# rngNoise  = MersenneTwister() # ruido
# rngInit   = MersenneTwister() # condiciones iniciales
# rngLinks  = MersenneTwister() # para red de Interaccion
#
# seedNoise = convert(Uint32,666)
# seedInit  = convert(Uint32,123)
# seedLinks = convert(Uint32,465)
#
# srand(rngNoise,seedNoise)
# srand(rngInit,seedInit)
# srand(rngLinks,seedLinks)

###========================================###

function rand_vec(w::Float64)
    # return Float64[-w + rand(rngInit)*2.0*w, -w + rand(rngInit)*2.0*w]
    return Float64[-w + rand()*2.0*w, -w + rand()*2.0*w]
end

###========================================###

function rand_num(w::Float64)
    # return -w + rand(rngInit)*2.0*w
    return -w + rand()*2.0*w
end

###========================================###

function set_noise(N::Int64)
    # noise = Float64[rand(rngNoise) for i = 1:N]
    noise = Float64[rand() for i = 1:N]

    # broadcast!(x -> -1.0*pi + 2.0*x*pi, noise, noise)
    broadcast!(x -> -pi + 2.0*x*pi, noise, noise)
    return noise
end

###========================================###
### FUNCTIONS FOR CREATING INTERACTION NET ###
###========================================###

###========================================###

function checkId(id::Int64, n::Int64)
    for i in 0:n-1
        if in(id, i*n+1:n*(i+1))
            return i+1
        end
    end
end

###========================================####

function setLinks!(m::Int64, n::Int64, numLinks::Array{Int64, 1}, adj_vec::Array{Int64, 1})
    for i in 1:m
        while true
            link = rand(1:n*n)
            # link = int(rand(rngLinks)*n*n)+1
            if adj_vec[link] < 1
                adj_vec[link] = 1
                numLinks[checkId(link, n)] += 1
                break
            end
        end
    end
end

###========================================####

function reArrange(m::Int64, n::Int64, numLinks::Array{Int64, 1}, adj_vec::Array{Int64, 1})

    Nij = zeros(Int64,m+n)
    posNumLinks = zeros(Int64,n)

    k = 1
    for i in 1:n
        Nij[k] = numLinks[i]
        posNumLinks[i] = k
        k += numLinks[i] + 1
    end

    links = find(x -> x==1, adj_vec)

    k = 0
    for i in 1:n
        for j in 1:numLinks[i]
            part = links[j+k] % n
            if part != 0 # part - 1 -> [0,n-1] for GPU
                Nij[posNumLinks[i]+j] = part
            else
                Nij[posNumLinks[i]+j] = n
            end
        end
        k += numLinks[i]
    end

    return Nij, posNumLinks
end

###========================================####

function make_IN(n::Int64, m::Int64)
    adj_vec = reshape(2*eye(Int64, n), n*n)
    numLinks = zeros(Int64, n)

    setLinks!(m, n, numLinks, adj_vec)
    return reArrange(m, n, numLinks, adj_vec)
end

###========================================###
### FUNCTION FOR PARTICLES INICIALIZATION  ###
###========================================###

function init_pos_vel(n::Int64, w::Float64, v0::Float64)

    vels = Float64[]

    for i in 1:n
        vi = rand_vec(1.0)
        scale!(vi, v0/norm(vi))
        vels = vcat(vels, vi)
    end

    return Float64[rand_num(w) for i in 1:2*n], vels
end

###========================================###
### SYSTEM EVOLUTION  ###
###========================================###

function calc_Rij(flock::Flock)

    nn = flock.n

    Rij = zeros(Float64,nn,nn)

    for i in 1:2:(2nn-2), j in (i+2):2:2nn
        Rij[div(i+1,2),div(j+1,2)] = Rij[div(j+1,2),div(i+1,2)] = (flock.pos[i] - flock.pos[j])^2 + (flock.pos[i+1] - flock.pos[j+1])^2
    end

    return Rij
end

###========================================####

function loc_vel(r0::Float64, Rij::Array{Float64,2}, flock::Flock)

    r02 = r0^2

    for part in 1:flock.n

        vx = flock.vels[2part-1]
        vy = flock.vels[2part]

        # loc_neigh = [find(x -> x <= r02 && x > 0.0, Rij[part,:])]
        loc_neigh = collect(find(x -> x <= r02 && x > 0.0, Rij[part,:]))

        k = 1.0 + Float64(length(loc_neigh))

        for j in 1:length(loc_neigh)
            vx += flock.vels[2j-1]
            vy += flock.vels[2j]
        end

        flock.VR[2part-1] = vx/k
        flock.VR[2part]   = vy/k
    end
end

###========================================####

function non_loc_vel(flock::Flock)

    for part in 1:flock.n

        ki = flock.Nij[flock.poski[part]] # links per particle
        init = flock.poski[part] # position to start reading

        vx = 0.0
        vy = 0.0

        if ki > 0

            for j in 1:ki
                vx += flock.vels[2 * flock.Nij[ init + j] - 1 ];
                vy += flock.vels[2 * flock.Nij[ init + j] ];
            end

            flock.VN[2part-1] = vx/float64(ki)
            flock.VN[2part]   = vy/float64(ki)
        else
            flock.VN[2part-1] = 0.0
            flock.VN[2part]   = 0.0
        end

    end
end

###========================================####

function rot_move(flock::Flock, noise::Array{Float64,1}, η::Float64, ω::Float64, w::Float64)

    for id in 1:flock.n

        prop_angle = atan2(flock.vels[2id], flock.vels[2id-1])

        i_vx = flock.VR[2id-1]
        i_vy = flock.VR[2id]

        if i_vx != 0.0 || i_vy != 0.0
            loc_angle = atan2(i_vy, i_vx) - prop_angle;
        else
            loc_angle = 0.0;
        end

        i_vx = flock.VN[2id-1]
        i_vy = flock.VN[2id]

        if i_vx != 0.0 || i_vy != 0.0
            non_loc_angle = atan2(i_vy, i_vx) - prop_angle;
        else
            non_loc_angle = 0.0;
        end

        total_angle = ω * loc_angle + (1 - ω) * non_loc_angle + η * noise[id];

        c = cos(total_angle)
        s = sin(total_angle)

        vx = flock.vels[2id-1]*c - flock.vels[2id]*s;
        vy = flock.vels[2id-1]*s + flock.vels[2id]*c;

        flock.vels[2id-1] = vx;
        flock.vels[2id]   = vy;

        flock.pos[2id-1] += vx;
        flock.pos[2id]   += vy;

        periodic_bounds(flock.pos[2id-1], flock.pos[2id], w)

    end
end
###========================================####

function evol(flock::Flock, r0::Float64, η::Float64, ω::Float64, w::Float64)
    loc_vel(r0, calc_Rij(flock), flock)
    non_loc_vel(flock)
    rot_move(flock, set_noise(flock.n), η, ω, w)
end
###========================================####

function periodic_bounds(posX::Float64, posY::Float64, w::Float64)

    sc = 0.5

    # for id in 1:flock.n

        if posX > w * sc
            posX -= w

        elseif posX <= -w * sc
            posX += w

        elseif posY > w * sc
            posY -= w

        elseif posY <= -w * sc
            posY += w

        end

    # end

end
###========================================####

end
