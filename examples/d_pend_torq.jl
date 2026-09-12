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

mot1 = PositionMotor2D(jnt2, 3.0, 0.0)

set_position_on_second_body!(jnt2, SA[-0.5, 0])

system = MBSystem2D()

add!(system, bd1)
add!(system, bd2)

add!(system, jnt1)
add!(system, jnt2)
add!(system, mot1)

if (!assemble!(system))
    println("Assembling failed!")
end

func = system.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)

bd1_x_ind, bd1_y_ind, _ = get_body_position_dofs(system, bd1)
bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(system, bd2)
initial = zeros(number_of_dofs(system))
initial[bd2_x_ind] = 0.5
system.prestep = (state) -> begin
    t = state[end]
    if t < 5
        settarget!(mot1,  π*t, 3.0)
    else 
        settarget!(mot1,  π*(10 - t), 3.0)
    end
end

time_start = 0
time_end = 10
time_step = 1

time_span = 0:0.001:10

sol2 = simulate(system, initial, time_span)
animate(system, sol2, time_span, "out/d_pend_torq.mp4"; framerate = floor(Int64, 1.0 / step(time_span)), limits = (-3,3, -3, 3))
##
jnt1_lbd = Flexia.get_lms(system, jnt1)
jnt2_lbd = Flexia.get_lms(system, jnt2)
mot1_lbd = Flexia.get_lms(system, mot1)



λ = sol2[[jnt1_lbd; jnt2_lbd; mot1_lbd], :]
inds = [3]
f, ax = series(time_span, λ[inds, :])