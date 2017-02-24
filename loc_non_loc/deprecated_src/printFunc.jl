## =========================== ##

# Martin Zumaya Hernandez
# Enero 2015
# Implementacion Flocks en Julia
# Funciones para salida de Datos

## =========================== ##

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

#Escribe Adjacencias de la Red de Interaccion
function PrintIntNet(path::ASCIIString,parts::Array{Bird,1})
    intNet = open("$path/intNet.txt","w")

    for i in 1:size(parts,1)
        write(intNet,repr(parts[i].inputs)[2:end-1]*"\n")
    end

    close(intNet)
end

## =========================== ## ## =========================== ##

#Escribe Archivo de Parametros
function PrintParams(path::ASCIIString)

  d = open("$path/params.txt","w")

  write(d,"particulas           = $N\n")
  write(d,"densidad             = $(params["p"])\n")
  write(d,"radio                = $r0\n");
  write(d,"conectividad         = $(params["k"])\n");
  write(d,"intensidad de ruido  = $(params["eta"])\n");
  write(d,"peso relativo        = $(params["w"])\n");
  write(d,"regimen de velocidad = $l\n");
  write(d,"iteraciones          = $(params["T"])\n");
  write(d,"frec muestreo        = $(params["s"])\n");
  write(d,"v0                   = $v0\n");
  write(d,"part Pert            = $partPert\n");
  write(d,"tau                  = $tau\n");
  write(d,"turnRate             = $turnRate\n");

  close(d)

end

## =========================== ## ## =========================== ##

function MakeDir(path::ASCIIString)
  try
    run(`mkdir $path`)
    run(`mkdir $path/dists`)
  catch y
    println(typeof(y))
  end
end

## =========================== ## ## =========================== ##

# ==================================== Funciones END ============================================
