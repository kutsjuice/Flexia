using Pkg; Pkg.activate(@__DIR__)
using Flexia
using CairoMakie
using ForwardDiff
using LinearAlgebra
using Polynomials
using StaticArrays

# ===========================================================================
# Drift of a hinge: revolute joint with a torsional spring and a prescribed
# initial velocity. The moving link orbits the hinge, and the example fits a
# cubic polynomial through the measured motion to locate the instantaneous
# centre of rotation and to compare it with the hinge point.
#
# Run from the package root (writes out/hinge_no_drift.mp4):
#
#     julia --project=examples examples/hinge_with_drift.jl

link_len = 0.15
link_width = 0.1
notch_width = 0.02
notch_depth = (link_width - notch_width)/2
contour = [
    Point2d(0.0, 0.0), 
    Point2d(0.0, link_width/2),
    Point2d(link_len - notch_width - notch_depth, link_width/2),
    Point2d(link_len - notch_width,  notch_width/2), 
    Point2d(link_len, notch_width/2), 
    Point2d(link_len, -notch_width/2), 
    Point2d(link_len-notch_width, -notch_width/2),
    Point2d(link_len - notch_width - notch_depth, -link_width/2),
    Point2d(0.0, -link_width/2),
    Point2d(0.0, 0.0), 
    ]

# simulation consist of two bodies: upper (fixed) and bottom (moving)
fixed = Body2D(10, 1.0; visual = MeshVisual(points = contour)) 
moving = Body2D(10.0, 1.0; visual = MeshVisual(points = contour))


# Create system
sys = MBSystem2D()
add!(sys, fixed)
add!(sys, moving)

ground_joint = FixedJoint(fixed, visual = FixedVisual(hatch_length = link_width, attachment = Attachment(angle=-π/2)))
setposition!(ground_joint, SA[0.0, 0.0])
setrotation!(ground_joint, -pi/2)


hinge = HingeJoint(fixed, moving)
set_position_on_first_body!(hinge, SA[link_len, 0.0])
set_position_on_second_body!(hinge, SA[link_len, 0.0])

stiffness = 0.5
damping = 0.0
spring = TorsionalSpring(hinge, stiffness, deg2rad(-180), damping, 0.03)

add!(sys, ground_joint)
add!(sys, hinge)
add!(sys, spring)


eop_x = FramePositionSensor(body=moving, position=SA_F64[0, 0.0 ], axis=:x)
eop_y = FramePositionSensor(body=moving, position=SA_F64[0, 0.0 ], axis=:y)
eop_r = FrameRotationSensor(moving)

add!(sys, eop_x)
add!(sys, eop_y)
add!(sys, eop_r)


assemble!(sys)

initial_state = zeros(number_of_dofs(sys))

set_initial_position!(initial_state, sys, fixed, SA_F64[0, 0, -pi/2])
set_initial_position!(initial_state, sys, moving, SA_F64[0, -2*link_len, pi/2])
set_initial_velocity!(initial_state, sys, moving, SA_F64[0, 0, 0.1])

time_span = 0:0.001:10.0
sol = simulate(sys, initial_state, time_span)

fig, ax, p = systemfigure(sys, sol)

fps = round(Int, 1 / (5 * step(time_span)))
record(fig, "out/hinge_no_drift.mp4", 1:5:length(time_span);
       framerate = fps) do t
    p.time = t
end


meas = similar(sol, (3, length(time_span)))
buf = meas[:, 1]
for i in axes(meas, 2)
    sys.measure!(buf, sol[:, i] )
    meas[:, i] .= buf;
end


u = meas .- meas[:,1]
Rot(th) = begin
    s,c = sincos(th)
    return SA[c -s; s c]
end

# Locate the instantaneous centre of rotation: a cubic fit of the relative
# displacement in the relative angle, then the point that maps both frames.
px = fit(u[3,:], u[1,:], 3)
py = fit(u[3,:], u[2,:], 3)
θ = LinRange(extrema(u[3,:])..., 100)
dpx = derivative(px)
dpy = derivative(py)

p = similar(meas, 2, length(time_span))
for i in axes(p,2)
    col = u[:, i]

    if abs(col[3]) < 1e-5
        p[1, i] = -dpy(col[3])
        p[2, i] = dpx(col[3])
    else
        p[:, i] .= inv(I(2) - Rot(col[3])) * col[1:2] 
    end
end

# The centre of rotation should hold still at the hinge point: the sample-wise
# estimate and the cubic fit are drawn together so a drift of the fit shows up
# as the distance between the two curves.
icr(θ) = SA[(px - dpy)(θ); (py + dpy)(θ)]

drift_fig = Figure()
drift_ax = Axis(drift_fig[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y")
lines!(drift_ax, p[1, :], p[2, :], label = "sample-wise estimate")
lines!(drift_ax, icr.(θ), label = "cubic fit")
axislegend(drift_ax)
drift_path = joinpath(@__DIR__, "..", "out", "hinge_with_drift.png")
save(drift_path, drift_fig)
println("saved ", drift_path)
