using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using StaticArrays

ground = Body2D(1e6, 1e6; length=0.01)  # Large mass/inertia to simulate fixed ground
const g = 9.81
crank_len = 1.0;
arm_len = 1.0;
H = 1.5;
L = 3.0
crank = Body2D(1, 1; length=crank_len)
arm = Body2D(1, 1; length=arm_len)

# crank.forces[2] = (x, t) -> -crank.mass * g
# arm.forces[2] = (x, t) -> -arm.mass * g


ground_joint = FixedJoint(ground)
setposition!(ground_joint, SA[0.0, 0.0])
setrotation!(ground_joint, 0)

hinge1 = HingeJoint(ground, crank)
set_position_on_first_body!(hinge1, SA[0.0, 0.0])
set_position_on_second_body!(hinge1, SA[-crank_len/2, 0])

hinge2 = HingeJoint(crank, arm)
set_position_on_first_body!(hinge2, SA[crank_len/2, 0])
set_position_on_second_body!(hinge2, SA[-arm_len/2, 0])

sys = MBSystem2D()
add!(sys, ground)
add!(sys, crank)
add!(sys, arm)

add!(sys, ground_joint)
add!(sys, hinge1)
add!(sys, hinge2)

assemble!(sys)
println("System assembled successfully")


init_state = zeros(number_of_dofs(sys))

set_initial_state_value!(init_state, sys, crank, SA_F64[crank_len/2, 0.0, 0.0])
set_initial_state_value!(init_state, sys, arm, SA_F64[crank_len, arm_len/2, pi/2])

draw_static(sys, init_state)

time_span = 0:0.01:10

# Solve
sol = simulate(sys, init_state, time_span)

sys.rhs(init_state)
sys.jacobian(init_state)

ngc = number_of_generalized_coordinates(sys)

p_state = zeros(ngc)
kin_constraints = zeros(ngc)

set_pos_state_value!(p_state, sys, crank, SA_F64[crank_len/2, 0.0, 0.0])
set_pos_state_value!(p_state, sys, arm, SA_F64[crank_len, arm_len/2, pi/2])
# sys.kinematic_residual!()

animate(sys, sol, time_span, "out/scara.mp4"; framerate= floor(Int, 1.0 /step(time_span)), limits=(-4, 5, -2, 2))


sys.kinematic_constrains!(kin_constraints, p_state)


