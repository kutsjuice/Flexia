using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff
using StaticArrays

const g = 9.81

ground = Body2D(1e6, 1e6; length=0.01)  # Large mass/inertia to simulate fixed ground
crank = Body2D(10.0, 1.0; length=1.0)
slider = Body2D(5.0, 0.1; length=0.25)

crank.forces[2] = (x, t) -> -crank.mass * g
slider.forces[2] = (x, t) -> -slider.mass * g

sys = MBSystem2D()
add!(sys, ground)
add!(sys, crank)
add!(sys, slider)

ground_joint = FixedJoint(ground)
setposition!(ground_joint, SA[0.0, 0.0])
setrotation!(ground_joint, 0.0)

slider_joint = SliderJoint(ground, slider)
set_position_on_first_body!(slider_joint, SA[0.0, 0.0])
set_position_on_second_body!(slider_joint, SA[0.0, 0.0])
set_direction_on_first_body!(slider_joint, SA[1.0, 0.0])
set_direction_on_second_body!(slider_joint, SA[1.0, 0.0])

hinge_joint = HingeJoint(slider, crank)
set_position_on_first_body!(hinge_joint, SA[0.0, 0.0])
set_position_on_second_body!(hinge_joint, SA[-crank.length/2, 0.0])

slider_actuator = PositionLinearActuator2D(slider_joint, 0.0)

# Add joints to system
add!(sys, ground_joint)
add!(sys, slider_joint)
add!(sys, hinge_joint)

add!(sys, slider_actuator)

initial_state = zeros(number_of_dofs(sys))

crank_θ = π/2
crank_x = crank.length/2 * cos(crank_θ)
crank_y = crank.length/2 * sin(crank_θ)
set_initial_position!(initial_state, sys, crank, SA_F64[crank_x, crank_y, crank_θ])
# set_initial_velocity!(initial_state, sys, crank, SA_F64[0, -5, 0])
#

measured = [
    get_body_position_dofs(sys, slider)[1], 
    get_body_velocity_dofs(sys, slider)[1], 
    get_body_position_dofs(sys, crank)[3], 
    get_body_velocity_dofs(sys, crank)[3]
    ]

sys.prestep = (state) -> begin
    y = state[mesured];
    t = state[end]
    
    settarget!(slider_actuator,  sin(t))
end




assemble!(sys)

##
# Flexia.update_targets!(sys)
# Flexia.get_targets(sys)[13:end]'
#
time_span = 0:0.0025:10
ratio = 1
sol = simulate(sys, initial_state, time_span)

animate(sys, sol, ratio*time_span, "out/inverted_pendulum.mp4"; framerate=1.0/(ratio * step(time_span)), limits=(-2, 2, -2, 2))

for_spy(mat) = reverse(mat', dims=2)