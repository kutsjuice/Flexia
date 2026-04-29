module Flexia
# using Pkg; Pkg.activate(".")
using StaticArrays
using Makie
using LinearAlgebra
# export MBSystem
# export Body2D


export Body2D
export FixedJoint, HingeJoint, SliderJoint
export PositionMotor2D
export MBSystem2D



export set_position_on_first_body!, set_position_on_second_body!, set_direction_on_first_body!, set_direction_on_second_body!, setposition!, setrotation!
export add!, assemble!, get_body_position_dofs, get_body_velocity_dofs, numberofdofs, last_body_dof, last_lm_dof, get_boundary_points
export simulate, cros!
export number_of_dofs
export animate
export get_mass_matrix

export test_func


test_func() = 1

abstract type AbstractBody2D end
abstract type AbstractConnector2D end
abstract type AbstractJoint2D <: AbstractConnector2D end
abstract type AbstractActuator2D <: AbstractConnector2D end
abstract type AbstractPositionActuator2D <: AbstractActuator2D end

include("system.jl")
include("solvers.jl")
include("bodies.jl")
include("connectors.jl")
include("joints.jl")
include("actuators.jl")
include("visualize.jl")

function getdofs(sys::MBSystem2D, body::Body2D)
    if (body.index == -1)
        return 0:0
    elseif (body.index == 1)
        return 1:sys.bodiesdofs[1]
    else
        return (sys.bodiesdofs[body.index-1]+1):sys.bodiesdofs[body.index]
    end
end


end
"""
# end # module Flexia

"""
