# module Flexia
using Pkg; Pkg.activate(".")
using StaticArrays
# using GLMakie
# export MBSystem
# export Body2D

abstract type AbstractBody2D end
abstract type AbstractJoint2D end

include("system.jl")

include("bodies.jl")
include("joints.jl")


function getdofs(sys::MBSystem2D, body::Body2D)
    if (body.index == -1)
        return 0:0
    elseif (body.index == 1)
        return 1:sys.bodiesdofs[1]
    else
        return (sys.bodiesdofs[body.index-1]+1):sys.bodiesdofs[body.index]
    end
end




# end # module Flexia
const g = 9.81

bd1 = Body2D(10, 1)
bd2 = Body2D(10, 1)
# bd2 = Body2D(10, 1)

bd1.forces[2] = (x) -> -bd1.mass * g
bd2.forces[2] = (x) -> -bd2.mass * g


jnt1 = FixedJoint(bd1)
jnt2 = HingeJoint(bd1, bd2)
set_position_on_second_body!(jnt2, SA[-0.5, 0])

sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)
add!(sys, jnt1)
add!(sys, jnt2)

if (!assemble!(sys))
    println("Assembling failed!")
end
# end # module Flexia
func = sys.rhs



using ForwardDiff

jacoby = (x) -> ForwardDiff.jacobian(func, x)

initial = zeros(number_of_dofs(sys))
func(initial)
jacoby(initial)