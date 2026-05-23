using Pkg; Pkg.activate("./examples")
using Flexia
using ForwardDiff
using GLMakie
using StaticArrays
# using CSV
# using DataFrames
const g = 9.81

m_ground = 1e6
I_ground = 1e6

m1 = 23e-3
m2 = 29e-3
I2 = 104726.798394e-9

m4 = 15e-3
I4 = 28187.202983e-9

K₁ = 0.2

Ground = Body2D(m_ground, I_ground, length = 0.15)
Ground.forces[2] = (x, t) -> -Ground.mass * g

bd2 = Body2D(m2+m1, I2, length = 0.075)
bd2.forces[2] = (x, t) -> -bd2.mass * g

slider = Body2D(m4, I4, length = 0.0255)
slider.forces[2] = (x, t) -> -slider.mass * g


crank = Body2D(m4, I4, length = 0.15)
crank.forces[2] = (x, t) -> -slider.mass * g



fix = FixedJoint(Ground)

jnt2 = HingeJoint(Ground, bd2)
set_position_on_first_body!(jnt2, SA[0., 0.])
set_position_on_second_body!(jnt2, SA[-bd2.length/2, 0.])

tcp = TorsionalSpring(Ground, bd2, K₁, deg2rad(0), 0.0, 0.03)

direcrion_SL = SA[1., 0.]

slider_joint = SliderJoint(bd2, slider)
set_position_on_first_body!(slider_joint, SA[0., 0.])
set_position_on_second_body!(slider_joint, SA[0., 0.])
set_direction_on_first_body!(slider_joint, direcrion_SL)
set_direction_on_second_body!(slider_joint, direcrion_SL)

hsp1 = HorizontalSpring(bd2, slider, 2., 0., 0.1, 0.1, 6)

hinge2 = HingeJoint(Ground, crank)
set_position_on_first_body!(hinge2, SA[0.1, -crank.length/2])
set_position_on_second_body!(hinge2, SA[-crank.length/2, 0.])


hinge3= HingeJoint(slider, crank)
set_position_on_first_body!(hinge3, SA[slider.length/2, 0])
set_position_on_second_body!(hinge3, SA[0, 0.])

sys = MBSystem2D()

add!(sys, Ground)
add!(sys, bd2)
add!(sys, slider)
add!(sys, crank)

add!(sys, fix)
add!(sys, jnt2)
add!(sys, tcp)
add!(sys, slider_joint)
add!(sys, hsp1)
add!(sys, hinge2)
add!(sys, hinge3)

if (!assemble!(sys))
    println("Assembling failed!")
end

func = sys.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)

initial = zeros(number_of_dofs(sys))

bd2_pos = get_body_position_dofs(sys, bd2)
initial[bd2_pos] .= [bd2.length/2, 0, 0]

slider_pos = get_body_position_dofs(sys, slider)
initial[slider_pos] .= [0.1-slider.length/2, 0, 0]

crank_pos = get_body_position_dofs(sys, crank)
θ = π/2-0.05
pos = [0.1, -crank.length/2] + [cos(θ) -sin(θ); sin(θ) cos(θ)]*[crank.length/2, 0] 
initial[crank_pos] .= [pos[1], pos[2], θ]
initial[end] = 0.0

jacoby(initial)

mass = get_mass_matrix(sys)
t_end = 2
time_span = LinRange(0,t_end, 1001)

n_steps = length(time_span)

sol1 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
cros!(sol1, initial, mass, func, jacoby, step(time_span))

animate(sys, sol1, time_span, "out/simple_slider.mp4"; framerate = floor(Int, 1/10/step(time_span)), limits = (-0.1, 0.25, -0.1, 0.1))

##
fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect())
for body in Flexia.bodies(sys)
    bar = Vector{Point2f}(undef, 2)
    bar2 = Vector{Point2f}(undef, 2)
    bar .= get_boundary_points(sys, body, initial)
    bar2 .= get_boundary_points(sys, body, sol1[:,2])
    lines!(ax, bar)
    lines!(ax, bar2, linestyle=:dash)
end
fig