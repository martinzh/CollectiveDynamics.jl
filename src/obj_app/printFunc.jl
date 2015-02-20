# ==================================== Funciones  ==============================================

# Funciones para imprimirs

## =========================== ## ## =========================== ##

#Escribe trayectorias
function PrintTrays(i::Int64,parts::Array{Bird,1})

    # write(trays,"$i\t")

    N = size(parts,1)

    line = ""

    for i in 1:N
        line *= replace(repr(parts[i].pos)[2:end-1],",",'\t')
        if i<N
            line *= "\t"
        end
    end

    line *= "\n"

    write(trays,line)

end

## =========================== ## ## =========================== ##

#Escribe velocidades
function PrintVels(i::Int64,parts::Array{Bird,1})
    
    # write(vels,"$i\t")

    N = size(parts,1)

    line = ""

    for i in 1:N
        line *= replace(repr(parts[i].vel)[2:end-1],",",'\t')
        if i<N
            line *= "\t"
        end
    end

    line *= "\n"

    write(vels,line)

end

## =========================== ## ## =========================== ##

#Escribe Matriz Distancias
function PrintDist(i::Int64,Dist::Array{Float64,2})
    writedlm("$path/dists/$i.txt",Dist,'\t')
end

## =========================== ## ## =========================== ##

#Escribe Archivo de Parametros
function PrintParams()
    # d = open("../$path/params.txt","w")
    d = open("$path/params.txt","w")

    write(d,"Particulas = $N\n")
    write(d,"densidad = $p\n")
    write(d,"radio = $r0\n");
    write(d,"conectividad = $k\n");
    write(d,"intensidad de ruido = $eta\n");
    write(d,"peso relativo = $w\n");
    write(d,"regimen de velocidad = $l\n");
    write(d,"iteraciones = $T\n");
    write(d,"frec muestreo = $step\n");
    write(d,"v0 = $v0\n");

    close(d)
end

## =========================== ## ## =========================== ##

function MakeDir()
    try
        # run(`mkdir ../$path`)
        # run(`mkdir ../$path/dists`)

        run(`mkdir $path`)
        run(`mkdir $path/dists`)

    catch y
        println(typeof(y))
    end
end

## =========================== ## ## =========================== ##

# ==================================== Funciones END ============================================