


f = ARGS[1]
N = ARGS[2]

# f = "NLOC_P_MET_3D"
# N = 512

data_path = joinpath(homedir(),"art_DATA",f,"DATA", "data_N_$(N)")

for (root, dirs, files) in walkdir(data_path)
    println("Files in $root")
    for file in files
        # println(joinpath(root, file)) # path to files 
        l = length(reinterpret(Float64, read(joinpath(root, file))))
        println(l)
    end
end
