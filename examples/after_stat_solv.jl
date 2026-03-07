
using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff
using JSON

using StaticArrays

function load_vector(filename::String)
    try
        data = open(filename, "r") do io
            JSON.parse(io)
        end
        println("✓ Vector loaded from $filename")
        return data
    catch e
        println("✗ Error loading: $e")
        return nothing
    end
end

const g = 9.81

const m1 = 10
const m2 = 10
const m3 = 10
const m4 = 10
const m5 = 10
const m6 = 10
const m7 = 10
const m8 = 10

bd1 = Body2D(m1, 500)
bd2 = Body2D(m2, 500)
bd3 = Body2D(m3, 500)
bd4 = Body2D(m4, 500)
bd5 = Body2D(m5, 500)
bd6 = Body2D(m6, 500)
bd7 = Body2D(m7, 500)
bd8 = Body2D(m8, 500)
# bd2 = Body2D(10, 1)

bd1.forces[2] = (x) -> -bd1.mass * g
bd2.forces[2] = (x) -> -bd2.mass * g
bd3.forces[2] = (x) -> -bd3.mass * g
bd4.forces[2] = (x) -> -bd4.mass * g
bd5.forces[2] = (x) -> -bd5.mass * g
bd6.forces[2] = (x) -> -bd6.mass * g
bd7.forces[2] = (x) -> -bd7.mass * g
bd8.forces[2] = (x) -> -bd8.mass * g

jnt1 = FixedJoint(bd1)

jnt2 = HingeJoint(bd1,bd2)
jnt3 = HingeJoint(bd2,bd3)
jnt4 = HingeJoint(bd3,bd4)
jnt5 = HingeJoint(bd4,bd5)
jnt6 = HingeJoint(bd5,bd6)
jnt7 = HingeJoint(bd6,bd7)
jnt8 = HingeJoint(bd7,bd8)

tcp2 = TorsionalSpring(bd2, bd3, 100000.,0.0, 0.)
tcp3 = TorsionalSpring(bd3, bd4, 100000.,0.0, 0.)
tcp4 = TorsionalSpring(bd4, bd5, 100000.,0.0, 0.)
tcp5 = TorsionalSpring(bd5, bd6, 100000.,0.0, 0.)
tcp6 = TorsionalSpring(bd6, bd7, 100000.,0.0, 0.)

set_position_on_second_body!(jnt2, SA[-1., 0])

set_position_on_first_body!(jnt3, SA[1., 0])
set_position_on_second_body!(jnt3, SA[-1., 0])

set_position_on_first_body!(jnt4, SA[1., 0])
set_position_on_second_body!(jnt4, SA[-1., 0])

set_position_on_first_body!(jnt5, SA[1., 0])
set_position_on_second_body!(jnt5, SA[-1., 0])

set_position_on_first_body!(jnt6, SA[1., 0])
set_position_on_second_body!(jnt6, SA[-1., 0])

set_position_on_first_body!(jnt7, SA[1., 0])
set_position_on_second_body!(jnt7, SA[-1., 0])

set_position_on_first_body!(jnt8, SA[1., 0])
set_position_on_second_body!(jnt8, SA[-1., 0])

sys2 = MBSystem2D()

add!(sys2, bd1)
add!(sys2, bd2)
add!(sys2, bd3)
add!(sys2, bd4)
add!(sys2, bd5)
add!(sys2, bd6)
add!(sys2, bd7)
add!(sys2, bd8)

add!(sys2, jnt1)
add!(sys2, jnt2)
add!(sys2, jnt3)
add!(sys2, jnt4)
add!(sys2, jnt5)
add!(sys2, jnt6)
add!(sys2, jnt7)
add!(sys2, jnt8)

add!(sys2, tcp2)
add!(sys2, tcp3)
add!(sys2, tcp4)
add!(sys2, tcp5)
add!(sys2, tcp6)

if (!assemble!(sys2))
    println("Assembling failed!")
end

func2 = sys2.rhs

jacoby2 = (x) -> ForwardDiff.jacobian(func2, x)

mass = zeros(number_of_dofs(sys2), number_of_dofs(sys2));
for i in 1:last_body_dof(sys2)
    mass[i, i] = 1
end

initial3 = zeros(number_of_dofs(sys2))
sol22 = zeros(number_of_dofs(sys2))
sol22 = load_vector("sol22.json")

initial3[bd2_x_ind] = sol22[bd2_x_ind]
initial3[bd2_y_ind] = sol22[bd2_y_ind]
initial3[bd2_t_ind] = sol22[bd2_t_ind]

initial3[bd3_x_ind] = sol22[bd3_x_ind]
initial3[bd3_y_ind] = sol22[bd3_x_ind]

initial3[bd4_x_ind] = sol22[bd4_x_ind]
initial3[bd4_y_ind] = sol22[bd4_y_ind]
initial3[bd4_t_ind] = sol22[bd4_t_ind]

initial3[bd5_x_ind] = sol22[bd5_x_ind]
initial3[bd5_y_ind] = sol22[bd5_y_ind]

initial3[bd6_x_ind] = sol22[bd6_x_ind]
initial3[bd6_y_ind] = sol22[bd6_y_ind]
initial3[bd6_t_ind] = sol22[bd6_t_ind]

initial3[bd7_x_ind] = sol22[bd7_x_ind]
initial3[bd7_y_ind] = sol22[bd7_y_ind]

initial3[bd8_x_ind] = sol22[bd8_x_ind]
initial3[bd8_y_ind] = sol22[bd8_y_ind]
initial3[bd8_t_ind] = sol22[bd8_t_ind]

func2(initial3)
jacoby2(initial3)

time_span = 0:0.0005:2

sol3 = Matrix{Float64}(undef, number_of_dofs(sys2), length(time_span))
cros!(sol3, initial3, mass, func2, jacoby2, step(time_span))
animate(sys, sol3, time_span, "atlet_diag_bar.mp4"; framerate = 60, limits = (-10,10, -10, 10))