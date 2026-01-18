using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff

using StaticArrays

const g = 9.81

bd1 = Body2D(10, 1)
bd2 = Body2D(10, 1)
# bd2 = Body2D(10, 1)

bd1.forces[2] = (x) -> -bd1.mass * g
bd2.forces[2] = (x) -> -bd2.mass * g


jnt1 = FixedJoint(bd1)
jnt2 = HingeJoint(bd1, bd2)

ac1 = PositionMotor2D(jnt2, π/6)

set_position_on_second_body!(jnt2, SA[-0.5, 0])

sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)

add!(sys, jnt1)
add!(sys, jnt2)

add!(sys, ac1)


sys.prestep = (t) -> begin
    ac1.target_angle += 0.01
end

if (!assemble!(sys))
    println("Assembling failed!")
end
# end # module Flexia
func = sys.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)

sys.jacobian = jacoby

bd1_x_ind, bd1_y_ind, _ = get_body_position_dofs(sys, bd1)
bd2_x_ind, bd2_y_ind, _ = get_body_position_dofs(sys, bd2)

initial = zeros(number_of_dofs(sys))
initial[bd2_x_ind] = 0.5
func(initial)
jacoby(initial)



mass = zeros(number_of_dofs(sys), number_of_dofs(sys));
for i in 1:last_body_dof(sys)
    mass[i, i] = 1
end
time_span = 0:0.01:10
sol = simulate(sys, initial, time_span)

animate(sys, sol, time_span, "rotor.mp4"; framerate = 30, limits = (-5,5, -5, 5))
