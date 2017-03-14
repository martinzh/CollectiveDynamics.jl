## =========================== ## ## =========================== ##
## 	   Package of functions of collective motion models			 ##
## 	   and its statistics                           			 ##
##	   Martin Zumaya Hernandez 						             ##
##     10 / 03 / 2017									         ##
## =========================== ## ## =========================== ##

module CollectiveDynamics

export InertialParameters, LocNonLocParameters, set_up_intertial_system!, set_up_intertial_nonLoc_system!, full_time_evolution_intertial_system, full_time_evolution_nonLocal_inertial_system, set_output_data_structure_inertial, set_output_data_structure_inertial_nonLocal, set_up_loc_nonLoc_system_2D!, set_up_loc_nonLoc_system_3D!, full_time_evolution_2D, full_time_evolution_3D, set_output_data_structure_2D, set_output_data_structure_3D

##  LOCAL + NON_LOCAL INTERACTIONS MODEL
include("local_nonLocal.jl")

##  INTERIAL TURN FLOCK MODEL
include("turn_flock.jl")

end
