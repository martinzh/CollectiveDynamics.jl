using CollectiveDynamics, HDF5, JLD, ArgParse

# ==================================== Parametros START ==========================================

s = ArgParseSettings()

@add_arg_table s begin
    "--steps"
        help     = "Total iterations"
        arg_type = Int64
        default  = 5000
    "--frec"
        help     = "Sampling Frequency"
        arg_type = Int64
        default  = 100
    "-n"
        help     = "Particles"
        arg_type = Int64
        default  = 2000
    "-m"
        help     = "Interaction Network's Connectivity"
        arg_type = Float64
        default  = 0.0
    "--eta"
        help     = "Noise Intensity"
        arg_type = Float64
        default  = 0.1
    # "-r"
    #     help    = "replica system"
    #     action  = :store_true
end

params = parse_args(s)

const dt = 1.0
const v0 = 1.0
const l  = 0.1 # Regimen de Velocidad
const ω  = 0.5 # peso relativo interacciones

r0 = (v0 * dt) / l

const n = params["n"] # Numero de particulas (entero)

params["r0"] = r0
params["v0"] = v0

const η = params["eta"]
# const m = int(params["m"] * n * (n-1))
const m = round(Int, params["m"] * n * (n-1))

###========================================###
##  CONIDICIONES INICIALES

ρ = 5*1e-2
params["p"] = ρ

const w = sqrt(n/ρ)

flock = Flock(n,m)

flock.pos, flock.vels  = init_pos_vel(n, w, 1.0)
flock.Nij, flock.poski = make_IN(n, m)

###========================================###

posArray  = Float64[]
velsArray = Float64[]

posArray  = vcat(posArray, flock.pos)
velsArray = vcat(velsArray, flock.vels)

path = "../DATA/data_n$(params["n"])_eta$(params["eta"])"

make_dir(path)

simData = jldopen("$path/m_$(params["m"]).jld","w")

println("End Setup")

###========================================###

for i in 1:params["steps"]

    evol(flock, r0, η, ω, w)

    if i%params["frec"] == 0

        println("t:",i)

        posArray  = vcat(posArray,flock.pos)
        velsArray = vcat(velsArray,flock.vels)
    end

end

###========================================###

posArray  = reshape(posArray,  2n, div(params["steps"], params["frec"])+1)
velsArray = reshape(velsArray, 2n, div(params["steps"], params["frec"])+1)

simData["pos"]    = posArray'
simData["vels"]   = velsArray'
simData["Nij"]    = flock.Nij
simData["params"] = params
simData["posKi"]  = flock.poski

close(simData)

println("DONE")
