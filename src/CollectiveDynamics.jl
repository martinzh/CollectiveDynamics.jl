## =========================== ## ## =========================== ##
## 	   Package of functions of collective motion models			 ##
## 	   and its statistics                           			 ##
##	   Martin Zumaya Hernandez 						             ##
##     10 / 03 / 2017									         ##
## =========================== ## ## =========================== ##

module CollectiveDynamics

export InertialParameters, LocNonLocParameters, set_up_intertial_system!, evolve_intertial_system, evolve_nonLocal_inertial_system, set_output_data_structure_inertial, set_output_data_structure_inertial_nonLocal, set_up_loc_nonLoc_system_2D!, set_up_loc_nonLoc_system_3D!, evolve_2D!, evolve_3D!, set_output_data_structure_2D, set_output_data_structure_3D

##  LOCAL + NON_LOCAL INTERACTIONS MODEL
include("local_nonLocal.jl")

##  INTERIAL TURN FLOCK MODEL
include("turn_flock.jl")

end
