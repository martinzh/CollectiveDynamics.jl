
# include("StatsLib.jl")
include("../src/StatsJulia/StatsLib.jl")

# path = "/home/martin/DATOS_SIMS/DATOS/data_f0.0"

path = "/Users/martinzh/DATOS_SIMS/DatJul/data_"

# println(path)

# Estructura parametros:
	# Particulas
	# densidad
	# radio
	# f
	# ruido geometrico
	# ruido topologico
	# peso geometrico
	# regimen de velocidad
	# iteraciones
	# paso

# obtiene parametros de archivo
params = Tmp.GetParams(path)
# print(params)

N        = int(params[0]) #Numero de particulas
ro       = params[1] #Densidad
r_0      = params[2] #Radio interaccion (vel_reg)
f        = params[3] #Fraccion de N en largo alcance
noise_sh = params[4] #Ruido corto
noise_lg = params[5] #Ruido largo
rel_weig = params[6] #Peso Relativo vecindades
reg_vel  = params[7] #Regimen de velocidad
t_f      = params[8] #iteraciones totales
step     = params[9] #frecuencia de muestreo

# step    = 250 #frecuencia de muestreo (Se obtiene del archivo de params)
num_bin = 200

bin_vec = arange(1,num_bin+1)

tiempo = arange(1,t_f+step,step)

# print(tiempo)

it_tot = t_f/step #muestras totales

# Tmp.GetParams(path)
