

f = ARGS[1]
N = parse(Int, ARGS[2])

# f = "NLOC_P_MET_3D"
# N = 512

tr_ln = 3N * 369

data_path = joinpath(homedir(),"art_DATA",f,"DATA", "data_N_$(N)")

for (root, dirs, files) in walkdir(data_path)
    println("Files in $root")
    for file in files
        # println(joinpath(root, file)) # path to files 
        ln = length(reinterpret(Float64, read(joinpath(root, file))))
        println(ln, "\t", ln < tr_ln)
    end
end
