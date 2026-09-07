using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff
using LinearAlgebra

using StaticArrays

const g = 9.81

bd1 = Body2D(1, 1)
bd2 = Body2D(1, 1)

# bd2.forces[2] = (x, t) -> -bd2.mass * g

jnt1 = FixedJoint(bd1)
jnt2 = HingeJoint(bd1, bd2)

mot1 = PositionMotor2D(jnt2, 0.0)

set_position_on_second_body!(jnt2, SA[-0.5, 0])

sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)

add!(sys, jnt1)
add!(sys, jnt2)
add!(sys, mot1)

if (!assemble!(sys))
    println("Assembling failed!")
end

func = sys.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)

bd1_x_ind, bd1_y_ind, _ = get_body_position_dofs(sys, bd1)
bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(sys, bd2)
initial = zeros(number_of_dofs(sys))
initial[bd2_x_ind] = 0.5
sys.prestep = (state) -> begin
    t = state[end]
    settarget!(mot1,  π*t)
end

time_start = 0
time_end = 10
time_step = 1

time_span = 0:0.001:10

sol2 = simulate(sys, initial, time_span)
animate(sys, sol2, time_span, "out/d_pend_torq.mp4"; framerate = floor(Int64, 1.0 / step(time_span)), limits = (-3,3, -3, 3))
##
jnt1_lbd = Flexia.get_lms(sys, jnt1)
jnt2_lbd = Flexia.get_lms(sys, jnt2)
mot1_lbd = Flexia.get_lms(sys, mot1)



λ = sol2[[jnt1_lbd; jnt2_lbd; mot1_lbd], :]
inds = [3]
f, ax = series(time_span, λ[inds, :])