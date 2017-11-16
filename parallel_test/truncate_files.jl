
include("par_mod.jl")

f = ARGS[1]
N = parse(Int, ARGS[2])

# f = "NLOC_P_MET_3D"
# N = 512

tr_ln = 3N * 369

data_path = joinpath(homedir(),"art_DATA",f,"DATA", "data_N_$(N)")

for (root, dirs, files) in walkdir(data_path)
    # println("Files in $root")
    for file in files

        data = reinterpret(Float64, read(joinpath(root, file)))
        ln = length(data)

        println(file,"\t", ln, "\t", tr_ln, "\t", ln < tr_ln)
        # if ln < tr_ln == false
        #     println(file,"\t", ln)
        #     write(joinpath(root,"cp_"*file), data[1:tr_ln])
        #     # rm(joinpath(root, file))
        #     # mv(joinpath(root,"cp_"*file),joinpath(root, file))
        # # else
        # #     rm(joinpath(root,file))
        # end

    end
end

println("Done")
