using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff
using LinearAlgebra

using StaticArrays

const g = 9.81

bd1 = Body2D(1, 1)
bd2 = Body2D(1, 1)

bd2.forces[2] = (x, t) -> -bd2.mass * g

jnt1 = FixedJoint(bd1)
jnt2 = HingeJoint(bd1, bd2)

set_position_on_second_body!(jnt2, SA[-0.5, 0])


eop_x = FramePositionSensor(body=bd2, position=SA_F64[0.5, 0.0 ], axis=:x)
eop_y = FramePositionSensor(body=bd2, position=SA_F64[0.5, 0.0 ], axis=:y)

eogv_x = FrameGlobalVelocitySensor(body=bd2, position=SA_F64[0.5, 0.0 ], axis=:x)
eogv_y = FrameGlobalVelocitySensor(body=bd2, position=SA_F64[0.5, 0.0 ], axis=:y)

eolv_x = FrameLocalVelocitySensor(body=bd2, position=SA_F64[0.5, 0.0 ], axis=:x)
eolv_y = FrameLocalVelocitySensor(body=bd2, position=SA_F64[0.5, 0.0 ], axis=:y)

eola_x = FrameLocalAccelerationSensor(body=bd2, position=SA_F64[0.5, 0.0 ], axis=:x)
eola_y = FrameLocalAccelerationSensor(body=bd2, position=SA_F64[0.5, 0.0 ], axis=:y)


sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)

add!(sys, jnt1)
add!(sys, jnt2)

add!(sys, eop_x)
add!(sys, eop_y)
add!(sys, eogv_x)
add!(sys, eogv_y)
add!(sys, eolv_x)
add!(sys, eolv_y)

add!(sys, eola_x)
add!(sys, eola_y)

if (!assemble!(sys))
    println("Assembling failed!")
end

func = sys.rhs

jacoby = sys.jacobian

bd1_x_ind, bd1_y_ind, _ = get_body_position_dofs(sys, bd1)
bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(sys, bd2)

initial = zeros(number_of_dofs(sys))
initial[bd2_x_ind] = 0.5
func(initial)
jacoby(initial)

mass = get_mass_matrix(sys)

time_start = 0
time_end = 10
time_step = 1

time_span = 0:0.001:10

sol2 = simulate(sys, initial, time_span)

meas = similar(sol2, (8, length(time_span)))
buf = meas[:, 1]
for i in axes(meas, 2)
    sys.measure!(buf, sol2[:, i] )
    meas[:, i] .= buf;
end

animate(sys, sol2, time_span, "out/simple_bar.mp4"; framerate = floor(Int64, 1.0 / step(time_span)), limits = (-2,2, -2, 2))

##
fig = Figure()
ax = Axis(fig[1, 1])

lines!(ax, time_span, meas[2,:])
lines!(ax, time_span, meas[4,:])
lines!(ax, time_span, meas[6,:])
lines!(ax, time_span, meas[8,:])
# lines!(ax, time_span, sol2[bd2_y_ind,:])
fig


