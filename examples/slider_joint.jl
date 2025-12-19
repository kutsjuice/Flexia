using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff
using StaticArrays

const g = 9.81

# Create two bodies: fixed base and sliding body
fixed_body = Body2D(100, 10)  # Heavy fixed body
sliding_body = Body2D(10, 1)  # Lighter sliding body

# Add gravity to the sliding body
sliding_body.forces[2] = (x) -> -sliding_body.mass * g

# Create joints
fixed_joint = FixedJoint(fixed_body)
slider_joint = SliderJoint(fixed_body, sliding_body)

# Set up the slider joint
# Position on fixed body (origin)
set_position_on_first_body!(slider_joint, SA[0.0, 0.0])
# Position on sliding body (center)
set_position_on_second_body!(slider_joint, SA[0.0, 0.0])

# Set the sliding direction: 20 degrees from horizontal
# Direction vector: (cos(θ), sin(θ)) where θ = 20° = π/9 radians
angle = π/9  # 20 degrees
direction = SA[cos(angle), sin(angle)]
set_direction_on_first_body!(slider_joint, direction)
set_direction_on_second_body!(slider_joint, direction)

# Create the system
sys = MBSystem2D()

# Add bodies
add!(sys, fixed_body)
add!(sys, sliding_body)

# Add joints
add!(sys, fixed_joint)
add!(sys, slider_joint)

# Assemble the system
if (!assemble!(sys))
    println("Assembling failed!")
end

# Set up the simulation
func = sys.rhs
jacoby = (x) -> ForwardDiff.jacobian(func, x)

# Get initial positions
fixed_x, fixed_y, _ = get_body_position_dofs(sys, fixed_body)
sliding_x, sliding_y, _ = get_body_position_dofs(sys, sliding_body)

# Initial state: all zeros (bodies start at origin)
initial = zeros(number_of_dofs(sys))

# Set initial position of sliding body slightly above the fixed body
initial[sliding_x] = 0.5  # Start 0.5 units to the right
initial[sliding_y] = 1.0  # Start 1 unit up

# Mass matrix
mass = zeros(number_of_dofs(sys), number_of_dofs(sys))
for i in 1:last_body_dof(sys)
    mass[i, i] = 1
end

# Time span
time_span = 0:0.01:5

# Solve
sol = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
cros!(sol, initial, mass, func, jacoby, step(time_span))

# Animate
animate(sys, sol, time_span, "slider_joint_animation.mp4"; framerate=30, limits=(-2, 2, -2, 2))
