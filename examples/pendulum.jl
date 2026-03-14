using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff

using StaticArrays

const g = 9.81

bd1 = Body2D(10, 500)
bd2 = Body2D(10, 500)
bd3 = Body2D(10, 500)
bd4 = Body2D(10, 500)
bd5 = Body2D(10, 500)
bd6 = Body2D(10, 500)
# bd2 = Body2D(10, 1)

bd1.forces[2] = (x) -> -bd1.mass * g
bd2.forces[2] = (x) -> -bd2.mass * g
bd3.forces[2] = (x) -> -bd3.mass * g
bd4.forces[2] = (x) -> -bd4.mass * g
bd5.forces[2] = (x) -> -bd5.mass * g
bd6.forces[2] = (x) -> -bd6.mass * g

jnt1 = FixedJoint(bd1)

jnt2 = HingeJoint(bd1, bd2)
jnt3 = HingeJoint(bd2, bd3)
jnt4 = HingeJoint(bd3, bd4)
jnt5 = HingeJoint(bd4, bd5)
jnt6 = HingeJoint(bd5, bd6)

jnt7 = FixedJoint(bd6)

tcp1 = TorsionalSpring(bd2, bd3, 1100.,0.0, 100.)
tcp2 = TorsionalSpring(bd3, bd4, 1100.,0.0, 100.)
tcp3 = TorsionalSpring(bd4, bd5, 1100.,0.0, 100.)

# traj = TrajectoryJoint(bd3, circular_trajectory([0.4, 1.1], 1.1, 1), 0.0, 10.0)

set_position_on_second_body!(jnt2, SA[-1., 0])

set_position_on_first_body!(jnt3, SA[1., 0])
set_position_on_second_body!(jnt3, SA[-1., 0])
set_position_on_first_body!(jnt4, SA[1., 0])
set_position_on_second_body!(jnt4, SA[-1., 0])
set_position_on_first_body!(jnt5, SA[1., 0])
set_position_on_second_body!(jnt5, SA[-1., 0])

set_position_on_first_body!(jnt6, SA[1., 0])

jnt7.pos = SA[2., 0.0]

sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)
add!(sys, bd3)
add!(sys, bd4)
add!(sys, bd5)
add!(sys, bd6)

add!(sys, jnt1)
add!(sys, jnt2)
add!(sys, jnt3)
add!(sys, jnt4)
add!(sys, jnt5)
add!(sys, jnt6)
add!(sys, jnt7)

add!(sys, tcp1)
add!(sys, tcp2)
add!(sys, tcp3)

# add!(sys, traj)
#mb1 = first_marker("mb1", bd1)

if (!assemble!(sys))
    println("Assembling failed!")
end
# end # module Flexia
func = sys.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)


bd1_x_ind, bd1_y_ind, _ = get_body_position_dofs(sys, bd1)
bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(sys, bd2)
bd3_x_ind, bd3_y_ind, bd3_t_ind = get_body_position_dofs(sys, bd3)
bd4_x_ind, bd4_y_ind, bd4_t_ind = get_body_position_dofs(sys, bd4)
bd5_x_ind, bd5_y_ind, bd5_t_ind = get_body_position_dofs(sys, bd5)
bd6_x_ind, bd6_y_ind, _ = get_body_position_dofs(sys, bd6)

bd2_Vx_ind, bd2_Vy_ind, bd2_Vt_ind = get_body_velocity_dofs(sys, bd2)
bd5_Vx_ind, bd5_Vy_ind, bd5_Vt_ind = get_body_velocity_dofs(sys, bd5)

bd3_Vx_ind, bd3_Vy_ind, bd3_Vt_ind = get_body_velocity_dofs(sys, bd3)


initial = zeros(number_of_dofs(sys))
# first configuretion
# initial[bd2_x_ind] = -0.171
# initial[bd2_y_ind] = 0.4698
# initial[bd2_t_ind] = 110*pi/180

# initial[bd3_x_ind] = 0.079
# initial[bd3_y_ind] = 1.2094
# initial[bd3_t_ind] = 32.39*pi/180

# initial[bd4_x_ind] = 0.921
# initial[bd4_y_ind] = 1.2094
# initial[bd4_t_ind] = -32.39*pi/180

# initial[bd5_x_ind] = 1.171
# initial[bd5_y_ind] = 0.4698
# initial[bd5_t_ind] = -110*pi/180

# initial[bd6_x_ind] = 1

 
 initial[bd2_Vt_ind] = -10
 initial[bd5_Vt_ind] = -10

#  zero config
# initial[bd2_x_ind] = 0.5
# initial[bd3_x_ind] = 1
# initial[bd3_y_ind] = 0.5
# initial[bd3_t_ind] = pi/2
# initial[bd4_x_ind] = 1
# initial[bd4_y_ind] = 0.5
# initial[bd4_t_ind] = -pi/2
# initial[bd5_x_ind] = 1.5
# initial[bd6_x_ind] = 2.0

#  second config

initial[bd2_y_ind] = 1
initial[bd2_t_ind] = pi/2

initial[bd3_x_ind] = 0.5
initial[bd3_y_ind] = 2 + sin(pi/3)
initial[bd3_t_ind] = pi/3

initial[bd4_x_ind] = 1.5
initial[bd4_y_ind] = 3
initial[bd4_t_ind] = - pi/3

initial[bd5_x_ind] = 2
initial[bd5_y_ind] = 1
initial[bd5_t_ind] = -pi/2

initial[bd6_x_ind] = 1

func(initial)
jacoby(initial)

mass = zeros(number_of_dofs(sys), number_of_dofs(sys));
for i in 1:last_body_dof(sys)
    mass[i, i] = 1
end

time_start = 0
time_end = 10
time_step = 1

time_span = 0:0.1:10

sol1 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
cros!(sol1, initial, mass, func, jacoby, step(time_span))

animate(sys, sol1, time_span, "time_animation16.mp4"; framerate = 60, limits = (-5,5, -5, 5))

sol2 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))

static_solver!( sol2 , initial, func, jacoby)

animate(sys, sol2, time_span, "static_5bar1.mp4"; framerate = 60, limits = (-5,5, -5, 5))

R1 = get_lms(sys, jnt1)
R2 = get_lms(sys, jnt2)
R3 = get_lms(sys, jnt3)
R4 = get_lms(sys, jnt4)
R5 = get_lms(sys, jnt5)
R6 = get_lms(sys, jnt6)
R7 = get_lms(sys, jnt7)

react_X1 = Float64[]
react_X2 = Float64[]
react_X3 = Float64[]
react_X4 = Float64[]
react_X5 = Float64[]
react_X6 = Float64[]
react_X7 = Float64[]

react_Y1 = Float64[]
react_Y2 = Float64[]
react_Y3 = Float64[]
react_Y4 = Float64[]
react_Y5 = Float64[]
react_Y6 = Float64[]
react_Y7 = Float64[]

react_M1 = Float64[]
react_M7 = Float64[]

react_TCP1 = Float64[]
react_TCP2 = Float64[]
react_TCP3 = Float64[]
for i in 1:time_end
    push!(react_X1, sol1[R1[1],i])
    push!(react_Y1, sol1[R1[2],i])
    push!(react_M1, sol1[R1[3],i])

    push!(react_X7, sol1[R7[1],i])
    push!(react_Y7, sol1[R7[2],i])
    push!(react_M7, sol1[R7[3],i])

    push!(react_X2, sol1[R2[1],i])
    push!(react_Y2, sol1[R2[2],i])

    push!(react_X3, sol1[R3[1],i])
    push!(react_Y3, sol1[R3[2],i])

    push!(react_X4, sol1[R4[1],i])
    push!(react_Y4, sol1[R4[2],i])

    push!(react_X5, sol1[R5[1],i])
    push!(react_Y5, sol1[R5[2],i])

    push!(react_X6, sol1[R6[1],i])
    push!(react_Y6, sol1[R6[2],i])

    push!(react_TCP1, get_spring_moment(tcp1, sol2[bd2_t_ind, i], sol2[bd3_t_ind, i]))
    push!(react_TCP2, get_spring_moment(tcp2, sol2[bd3_t_ind, i], sol2[bd4_t_ind, i]))
    push!(react_TCP3, get_spring_moment(tcp2, sol2[bd4_t_ind, i], sol2[bd5_t_ind, i]))
end

f1 = Figure()
f2 = Figure()
f3 = Figure()
time_span2 = 1:1:time_end

react_X1

ax1 = Axis(f1[1,1], xscale = identity, title="Reactions on Y",xlabel="Time, с.",ylabel="Force, Н.")
ax2 = Axis(f2[1,1], xscale = identity, title="Reactions on X",xlabel="Time, с.",ylabel="Force, Н.")
ax3 = Axis(f3[1,1], xscale = identity, title="Torques",xlabel="Time, с.",ylabel="Torque, Нxм.")

l1Y = lines!(ax1, time_span2, react_Y1, linestyle = :solid)
l2Y = lines!(ax1, time_span2, react_Y2, linestyle = :dot)
l3Y = lines!(ax1, time_span2, react_Y3, linestyle = :dash)
l4Y = lines!(ax1, time_span2, react_Y4, linestyle = :dashdot)
l5Y = lines!(ax1, time_span2, react_Y5, linestyle = :dashdotdot)
l6Y = lines!(ax1, time_span2, react_Y6, linestyle = :solid)
l7Y = lines!(ax1, time_span2, react_Y7, linestyle = :dash)

l1X = lines!(ax2, time_span2, react_X1, linestyle = :solid)
l2X = lines!(ax2, time_span2, react_X2, linestyle = :dot)
l3X = lines!(ax2, time_span2, react_X3, linestyle = :dash)
l4X = lines!(ax2, time_span2, react_X4, linestyle = :dashdot)
l5X = lines!(ax2, time_span2, react_X5, linestyle = :dashdotdot)
l6X = lines!(ax2, time_span2, react_X6, linestyle = :solid)
l7X = lines!(ax2, time_span2, react_X7, linestyle = :dash)

l1M = lines!(ax3, time_span2, react_M1, linestyle = :solid)
l7M = lines!(ax3, time_span2, react_M7, linestyle = :dot)
l1TCP = lines!(ax3, time_span2, react_TCP1, linestyle = :dash)
l2TCP = lines!(ax3, time_span2, react_TCP2, linestyle = :dashdot)
l3TCP = lines!(ax3, time_span2, react_TCP3, linestyle = :dashdotdot)

Legend(f1[1 , 2], [l1X,l2X,l3X, l4X, l5X, l6X, l7X], ["R1(Fix)", "R2(Jnt2)", "R3(Jnt3)", "R4(Jnt4)", "R5(Jnt5)", "R6(Jnt6)", "R7(Fix)"], framevisible = false, halign = :right, valign = :top)
Legend(f2[1 , 2], [l1Y,l2Y,l3Y, l4Y, l5Y, l6Y, l7Y], ["R1(Fix)", "R2(Jnt2)", "R3(Jnt3)", "R4(Jnt4)", "R5(Jnt5)", "R6(Jnt6)", "R7(Fix)"], framevisible = false, halign = :right, valign = :top)
Legend(f3[1 , 2], [l1M,l7M,l1TCP,l2TCP,l3TCP], ["M1(Fix)", "M7(fix)", "TCP1", "TCP2", "TCP3"], framevisible = false, halign = :left, valign = :top)

save("pendulum_reactions_Y.png", f1)
save("pendulum_reactions_X.png", f2)
save("pendulum_reactions_M.png", f3)