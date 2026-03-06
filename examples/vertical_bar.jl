using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff

using StaticArrays
using LinearAlgebra

const g = 9.81

const m1 = 10
const m2 = 10
const m3 = 10
const m4 = 10

bd1 = Body2D(m1, 100)
bd2 = Body2D(m2, 100)
bd3 = Body2D(m3, 100)
bd4 = Body2D(m4, 100)

# bd2 = Body2D(10, 1)

bd1.forces[2] = (x) -> -bd1.mass * g
bd2.forces[2] = (x) -> -bd2.mass * g
bd3.forces[2] = (x) -> -bd3.mass * g
bd4.forces[2] = (x) -> -bd3.mass * g

jnt1 = FixedJoint(bd1)
jnt2 = SpringY(bd1, bd2, 12., 0.0, 0.0)
jnt3 = SpringY(bd2, bd3, 12., 0.0, 0.0)
jnt4 = SpringY(bd3, bd4, 12., 0.0, 0.0)

set_rotation!(jnt1, -pi/2)

set_position_on_second_body!(jnt2, SA[-1., 0])

set_position_on_first_body!(jnt3, SA[1., 0])
set_position_on_second_body!(jnt3, SA[-1., 0])
set_position_on_first_body!(jnt4, SA[1., 0])
set_position_on_second_body!(jnt4, SA[-1., 0])

sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)
add!(sys, bd3)
add!(sys, bd4)

add!(sys, jnt1)
add!(sys, jnt2)
add!(sys, jnt3)
add!(sys, jnt4)

if (!assemble!(sys))
    println("Assembling failed!")
end
# end # module Flexia
func = sys.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)

bd1_x_ind, bd1_y_ind, bd1_t_ind = get_body_position_dofs(sys, bd1)
bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(sys, bd2)
bd3_x_ind, bd3_y_ind, bd3_t_ind = get_body_position_dofs(sys, bd3)
bd4_x_ind, bd4_y_ind, bd4_t_ind = get_body_position_dofs(sys, bd4)

initial = zeros(number_of_dofs(sys))
initial[bd1_t_ind] = -pi/2

initial[bd2_y_ind] = -2.
initial[bd2_t_ind] = -pi/2

initial[bd3_y_ind] = -3.
initial[bd3_t_ind] = -pi/2

initial[bd4_y_ind] = -4.
initial[bd4_t_ind] = -pi/2

func(initial)
jacoby(initial)
mass = zeros(number_of_dofs(sys), number_of_dofs(sys));
for i in 1:last_body_dof(sys)
    mass[i, i] = 1
end
time_span = 0:0.005:20

sol1 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
cros!(sol1, initial, mass, func, jacoby, step(time_span))

animate(sys, sol1, time_span, "vertical_bar.mp4"; framerate = 60, limits = (-5,5, -9, 1))
