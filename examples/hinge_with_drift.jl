using Pkg; Pkg.activate(@__DIR__)
using Flexia
using CairoMakie
using ForwardDiff
using LinearAlgebra
using Polynomials
using StaticArrays

# ===========================================================================
# Drift of a hinge under a static sweep: the joint angle is prescribed by a
# positional motor and the equilibrium at every pseudotime step is found by the
# static solver. The moving link swings ±25° around the relative angle of the
# initial pose, one full period over the run. The measured motion is fitted with
# a cubic and the resulting instantaneous centre of rotation is compared with
# the hinge point.
#
# Run from the package root (writes out/hinge_no_drift.mp4 and
# out/hinge_with_drift.png):
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

ground_joint = FixedJoint(fixed,
                          visual = FixedVisual(hatch_length = link_width,
                                               attachment = Attachment(angle = -π/2)))
setposition!(ground_joint, SA[0.0, 0.0])
setrotation!(ground_joint, -pi/2)


hinge = HingeJoint(fixed, moving)
set_position_on_first_body!(hinge, SA[link_len, 0.0])
set_position_on_second_body!(hinge, SA[link_len, 0.0])

# The motor prescribes the joint angle: a sine of 25° amplitude, one period over
# the whole run, around the relative angle of the initial pose.
angle0 = π
amplitude = deg2rad(25.0)
time_span = 0:0.001:10.0
ω = 2π / time_span[end]

angle_at(t) = angle0 + amplitude * sin(ω * t)
rate_at(t) = amplitude * ω * cos(ω * t)

motor = PositionMotor2D(hinge, angle_at(0.0), rate_at(0.0))

add!(sys, ground_joint)
add!(sys, hinge)
add!(sys, motor)


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

# Statics: the last slot of the state is the pseudotime, and the motor target is
# a function of it, so the target is refreshed before every residual evaluation.
# ForwardDiff sees that slot as a dual number, hence value(). `update_targets!` is
# the hook the integrator calls for us; the static path calls it explicitly.
function equilibrium(state)
    t = ForwardDiff.value(state[end])
    settarget!(motor, angle_at(t), rate_at(t))
    Flexia.update_targets!(sys)
    return sys.rhs(state)
end

equilibrium_jacobian(state) = ForwardDiff.jacobian(equilibrium, state)

# static_solver! puts the warm start into the first column and the load path into
# the rest, so the animation starts from the second column.
static = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span) + 1)
path = static_solver!(static, initial_state, equilibrium, equilibrium_jacobian;
                      times = time_span)
sol = static[:, 2:end]

converged = count(entry -> entry.status === :converged, path)
println("statics: $converged/$(length(path)) steps converged, ",
        "max residual = ", maximum(entry -> entry.residual, path))

Rot(th) = begin
    s,c = sincos(th)
    return SA[c -s; s c]
end


meas = similar(sol, (3, length(time_span)))
buf = meas[:, 1]
for i in axes(meas, 2)
    sys.measure!(buf, sol[:, i] )
    meas[:, i] .= buf;
end


u = meas .- meas[:,1]

# Locate the instantaneous centre of rotation: a cubic fit of the relative
# displacement in the relative angle. For a pure rotation u + Ωu' is exactly the
# pole, so the fitted ICR is (u_x - u_y', u_y + u_x').
ux = fit(u[3,:], u[1,:], 3)
uy = fit(u[3,:], u[2,:], 3)
θ = LinRange(extrema(u[3,:])..., 100)
dux = derivative(ux)
duy = derivative(uy)

icr(θ) = SA[(ux - duy)(θ); (uy + dux)(θ)]


# Hinge area only: the FRP (pole of the finite rotation from the initial pose)
# and the ICR (instantaneous centre of rotation) are recomputed on every frame
# from the current state. For an ideal hinge both stay on the hinge point.
dofs = get_body_position_dofs(sys, moving)
x0, y0, θ0 = meas[:,1]

frp_x = Observable(0.0); frp_y = Observable(-link_len)
icr_x = Observable(0.0); icr_y = Observable(-link_len)
##
fig, ax, p = systemfigure(sys, sol;
    limits = (-0.0025, 0.0025, -link_len - 0.0025, -link_len + 0.0025))
scatter!(ax, frp_x, frp_y, color = :transparent, strokecolor = :red, strokewidth = 3,
         markersize = 16, label = "ПКП",  marker=:circle)
scatter!(ax, icr_x, icr_y, color = :orange, markersize = 10, label = "МЦС",marker='X')
axislegend(ax; position = :lt, framevisible = false, labelsize = 12)
icr(0)
fps = round(Int, 1 / (20 * step(time_span)))
record(fig, "out/hinge_no_drift.mp4", 1:20:length(time_span);
       framerate = fps) do i
    p.time = i
    x, y, θ = sol[dofs, i]

    θr = θ - θ0
    icr_xy = icr(θr)
    icr_x[] = icr_xy[1]+x0
    icr_y[] = icr_xy[2]+y0
    d = inv(I(2) - Rot(θr)) * SA[x - x0, y - y0]
    frp_x[] = x0 + d[1]
    frp_y[] = y0 + d[2]

end
##
# Sample-wise estimate: the pole of the finite rotation of every frame, taken
# from the state itself, without any fit.
p = similar(meas, 2, length(time_span))
for i in axes(p,2)
    col = u[:, i]

    if abs(col[3]) < 1e-5
        p[1, i] = -duy(col[3])
        p[2, i] = dux(col[3])
    else
        p[:, i] .= inv(I(2) - Rot(col[3])) * col[1:2]
    end
end

# Both curves in the axes of the hinge point: an ideal hinge keeps them at (0, 0).
# The sample-wise estimate is ill-conditioned at θ → 0 (it amplifies the residual
# of the constraint by 1/θ), hence the scatter along the curves.
hinge_point = SA[0.0, link_len]     # hinge point relative to the sensor, initial pose
icr_pts = reduce(hcat, icr.(θ))

drift_fig = Figure()
drift_ax = Axis(drift_fig[1, 1], aspect = DataAspect(),
                xlabel = "Δx, mm", ylabel = "Δy, mm",
                title = "ICR relative to the hinge point")
sample_x = 1e3 .* (p[1, :] .- hinge_point[1])
sample_y = 1e3 .* (p[2, :] .- hinge_point[2])
fit_x = 1e3 .* (icr_pts[1, :] .- hinge_point[1])
fit_y = 1e3 .* (icr_pts[2, :] .- hinge_point[2])
lines!(drift_ax, sample_x, sample_y, label = "sample-wise estimate")
lines!(drift_ax, fit_x, fit_y, label = "cubic fit")
scatter!(drift_ax, [0.0], [0.0], color = :black, markersize = 8, label = "hinge point")
limits!(drift_ax, -0.05, 0.05, -0.05, 0.05)
axislegend(drift_ax)
drift_path = joinpath(@__DIR__, "..", "out", "hinge_with_drift.png")
save(drift_path, drift_fig)
println("saved ", drift_path)
