
include("par_mod.jl")

f = ARGS[1]
N = parse(Int, ARGS[2])

# f = "NLOC_P_MET_3D"
# N = 512

tr_ln = 3N * 369

data_path = joinpath(homedir(),"art_DATA",f,"DATA", "data_N_$(N)")

output_path = set_output_data_structure_lnl(f*"_TRUNC", N, κ, ω)

for (root, dirs, files) in walkdir(data_path)
    println("Files in $root")
    for file in files
        # println(joinpath(root, file)) # path to files 
        data = reinterpret(Float64, read(joinpath(root, file)))
        ln = length(data)
        # println(file, "\t", ln, "\t", ln < tr_ln)

        if ln < tr_ln == false
            println(file)
            write(joinpath(output_path,file), data[1:tr_ln])
        end
    end
end

println("Done")
