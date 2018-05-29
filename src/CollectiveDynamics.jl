## =========================== ## ## =========================== ##
## 	   Package of functions of collective motion models			 ##
##	   Martin Zumaya Hernandez 						             ##
##	   Release Version 0.1.0 						             ##
##     29 / 05 / 2018									         ##
## =========================== ## ## =========================== ##

module CollectiveDynamics

## SIMPLE VICSEK MODEL 2D
module SVM2D
export Box, Flock, set_output_data_structure_vsk, evolve_system
include("simple_vicsek_model_2D.jl")
end

## SIMPLE VICSEK MODEL 3D
module SVM3D
using Quaternions
export Box, Flock, set_output_data_structure_vsk, evolve_system
include("simple_vicsek_model_3D.jl")
end

## BEHAVIOURAL RULES MODEL
# module BehaviouralRules
# using Quaternions
# export
# include("behavioural_rules_model.jl")
# end
#
# ## INTERIAL SPIN MODEL MODEL
# module InertialSpin
# using Quaternions
# export
# include("inertial_spin_model.jl")
# end
#
# ## SHORT AND LONG RANGE INTERACTIONS MODEL
# module ShortLongRange
# using Quaternions
# export
# include("short_long_range_model.jl")
# end

end
