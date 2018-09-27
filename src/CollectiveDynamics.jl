## =========================== ## ## =========================== ##
## 	   Package of collective motion models			 ##
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
module BehaviouralRules
using Quaternions
export set_output_data_structure, evolve_system
include("behavioural_rules_model.jl")
end

## SHORT AND LONG RANGE INTERACTIONS MODEL
module ShortLongRange
using Distributions, Quaternions
export set_output_data_structure, evolve_metric_system, evolve_metric_system_2D, evolve_topological_system, evolve_topological_system_2D
include("short_long_range_model.jl")
end

## INTERIAL SPIN MODEL MODEL
module InertialSpin
using Quaternions
export InertialFlock, InertialExtFlock, InertialParameters, set_output_data_structure, set_output_data_structure_lr, evolve_system, evolve_extended_system
include("inertial_spin_model.jl")
end

## INTERIAL SPIN MODEL MODEL
module NetworkModel
export set_output_data_structure, set_directed_random_network, evolve_system_2D
include("network_model.jl")
end


end
