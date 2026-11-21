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

Rot(th) = begin
    s,c = sincos(th)
    return SA[c -s; s c]
end

# Hinge area only: the ПКП (pole of the finite rotation from the initial pose)
# and the МЦС (instantaneous centre of velocity) are recomputed on every frame
# from the current state. For an ideal hinge both stay on the hinge point.
dofs = get_body_position_dofs(sys, moving)
vdofs = get_body_velocity_dofs(sys, moving)
x0, y0, θ0 = sol[dofs, 1]

pkp_x = Observable(0.0); pkp_y = Observable(-link_len)
mcs_x = Observable(0.0); mcs_y = Observable(-link_len)

fig, ax, p = systemfigure(sys, sol;
    limits = (-0.025, 0.025, -link_len - 0.025, -link_len + 0.025))
scatter!(ax, pkp_x, pkp_y, color = :transparent, strokecolor = :red, strokewidth = 3,
         markersize = 16, label = "ПКП")
scatter!(ax, mcs_x, mcs_y, color = :blue, markersize = 7, label = "МЦС")
axislegend(ax; position = :lt, framevisible = false, labelsize = 12)

fps = round(Int, 1 / (5 * step(time_span)))
record(fig, "out/hinge_no_drift.mp4", 1:5:length(time_span);
       framerate = fps) do i
    p.time = i
    x, y, θ = sol[dofs, i]
    vx, vy, ω = sol[vdofs, i]
    mcs_x[] = x - vy / ω                     # точка звена с нулевой скоростью
    mcs_y[] = y + vx / ω
    θr = θ - θ0
    if abs(θr) < 1e-5                        # предел θ → 0: ПКП переходит в МЦС
        pkp_x[] = mcs_x[]
        pkp_y[] = mcs_y[]
    else
        d = inv(I(2) - Rot(θr)) * SA[x - x0, y - y0]
        pkp_x[] = x0 + d[1]
        pkp_y[] = y0 + d[2]
    end
end


meas = similar(sol, (3, length(time_span)))
buf = meas[:, 1]
for i in axes(meas, 2)
    sys.measure!(buf, sol[:, i] )
    meas[:, i] .= buf;
end


u = meas .- meas[:,1]

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
# as the distance between the two curves. For a pure rotation u + Ωu' is exactly
# the pole, so the fitted ICR is (u_x - u_y', u_y + u_x').
icr(θ) = SA[(px - dpy)(θ); (py + dpx)(θ)]

# Both curves in the axes of the hinge point: an ideal hinge keeps them at (0, 0).
# The sample-wise estimate is ill-conditioned at θ → 0 (it amplifies the constraint
# error by 1/θ), hence the ±0.02 mm scatter along x.
hinge = SA[0.0, link_len]        # hinge point relative to the sensor, initial pose
icr_pts = reduce(hcat, icr.(θ))

drift_fig = Figure()
drift_ax = Axis(drift_fig[1, 1], aspect = DataAspect(),
                xlabel = "Δx, mm", ylabel = "Δy, mm",
                title = "ICR relative to the hinge point")
lines!(drift_ax, 1e3 .* (p[1, :] .- hinge[1]), 1e3 .* (p[2, :] .- hinge[2]),
       label = "sample-wise estimate")
lines!(drift_ax, 1e3 .* (icr_pts[1, :] .- hinge[1]), 1e3 .* (icr_pts[2, :] .- hinge[2]),
       label = "cubic fit")
scatter!(drift_ax, [0.0], [0.0], color = :black, markersize = 8, label = "hinge point")
limits!(drift_ax, -0.05, 0.05, -0.05, 0.05)
axislegend(drift_ax)
drift_path = joinpath(@__DIR__, "..", "out", "hinge_with_drift.png")
save(drift_path, drift_fig)
println("saved ", drift_path)
