using Pkg; Pkg.activate(@__DIR__)
using Flexia
using CairoMakie
using ForwardDiff
using LinearAlgebra
using StaticArrays

# ===========================================================================
# Hinge with drift under a static sweep: the joint angle is prescribed by a
# positional motor and the equilibrium at every pseudotime step is found by the
# static solver. The axis of the hinge is not fixed — it slides along the analytic
# centroids, so the instantaneous centre of rotation (МЦС) is the common point of
# the fixed and the moving centroid and travels along the fixed one. The moving
# link swings ±25° around the relative angle of the initial pose, one full period
# over the run.
#
# The animation draws both centroids, the trajectory of the МЦС and the ПКП (the
# pole of the finite rotation from the initial pose, a different point); the
# nominal (undrifted) hinge point is marked for comparison.
#
# Run from the package root (writes out/hinge_with_drift.mp4 and
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

sys = MBSystem2D()
add!(sys, fixed)
add!(sys, moving)

ground_joint = FixedJoint(fixed,
                          visual = FixedVisual(hatch_length = link_width,
                                               attachment = Attachment(angle = -π/2)))
setposition!(ground_joint, SA[0.0, 0.0])
setrotation!(ground_joint, -pi/2)

# Centroids of the pair: c1 is the fixed one (body 1), c2 the moving one (body 2),
# both parameterized by the relative angle θ. θ is measured from the initial pose
# (rot_offset = π: at the start the moving link is turned by π against the fixed
# one), where the two curves touch at the nominal hinge point and roll without
# slipping, so the hinge slides exactly along the МЦС.
#
#     c1(θ) = [k(θ + sin θ),  k(1 − cos θ)]
#     c2(θ) = [k(θ + sin θ), −k(1 − cos θ)]
#
# The links stand vertically, so both equations are turned by +90° in the body
# frames: the МЦС of the fixed frame then runs along c1, mostly horizontally.
# The moving link is anti-parallel to the fixed one, which flips the sign of its
# centroid — otherwise the pair would not roll (the contact would not be the МЦС).
k = 0.02
C1(θ) = SA[k * (θ + sin(θ)),  k * (1 - cos(θ))]      # уравнения из ТЗ
C2(θ) = SA[k * (θ + sin(θ)), -k * (1 - cos(θ))]
rot90(v) = SA[-v[2], v[1]]
c1(θ) = rot90(C1(θ))       # в ЛСК тела 1
c2(θ) = rot90(-C2(θ))      # в ЛСК тела 2

angle0 = Float64(π)        # относительный угол исходной позы: отсчёт дрейфа
hinge = HingeJoint(fixed, moving; rot_offset = angle0, drift1 = c1, drift2 = c2)
set_position_on_first_body!(hinge, SA[link_len, 0.0])    # торец неподвижного звена
set_position_on_second_body!(hinge, SA[link_len, 0.0])   # торец подвижного звена

# The motor prescribes the joint angle: a sine of 25° amplitude, one period over
# the whole run, around the relative angle of the initial pose. The centroids see
# the same angle minus rot_offset.
amplitude = deg2rad(25.0)
time_span = 0:0.001:10.0
ω = 2π / time_span[end]

angle_at(t) = angle0 + amplitude * sin(ω * t)
rate_at(t) = amplitude * ω * cos(ω * t)

motor = PositionMotor2D(hinge, angle_at(0.0), rate_at(0.0))

add!(sys, ground_joint)
add!(sys, hinge)
add!(sys, motor)

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

##
# The fixed centroid in the world frame, from its analytic equation: in the world
# the +90° turn of the body frame cancels, so the МЦС runs along c1 of the ТЗ.
nominal = Point2f(0.0, -link_len)                 # номинальный шарнир в мире
θs = range(-amplitude, amplitude; length = 241)
centroid(θ) = Point2f(nominal[1] + C1(θ)[1], nominal[2] + C1(θ)[2])
fixed_centroid_pts = centroid.(θs)

q2 = get_body_position_dofs(sys, moving)

# The moving centroid in the world frame, for the current pose of the moving link.
# The curve is attached to the hinge point of that link, not to its origin.
h2 = SA[link_len, 0.0]
function moving_centroid_pts(state)
    x, y, φ = state[q2]
    c, s = cos(φ), sin(φ)
    return [Point2f(x + c * (h2[1] + c2(θ)[1]) - s * (h2[2] + c2(θ)[2]),
                    y + s * (h2[1] + c2(θ)[1]) + c * (h2[2] + c2(θ)[2])) for θ in θs]
end

# The marker of the МЦС of every frame comes from the drift-aware hinge geometry,
# so it must land on the common point of both centroids.
mcs_pts = [Flexia.get_hinge_point(sys, hinge, view(sol, :, i)) for i in axes(sol, 2)]
travel = 1e3 * (maximum(first.(mcs_pts)) - minimum(first.(mcs_pts)))
println("МЦС travels ", round(travel; digits = 1), " mm along the fixed centroid")

# ПКП — полюс конечного поворота подвижного звена от исходной позы. Это другая
# точка: при дрейфе она не совпадает с МЦС. При θ → 0 вырождается, предел равен
# мгновенному центру, то есть МЦС.
x2 = sol[q2[1], :]; y2 = sol[q2[2], :]; φ2 = sol[q2[3], :]
x0, y0, φ0 = x2[1], y2[1], φ2[1]

function frp_point(state)
    x, y, φ = state[q2]
    Δφ = φ - φ0
    abs(Δφ) < 1e-6 && return Flexia.get_hinge_point(sys, hinge, state)

    R = SA[cos(Δφ) -sin(Δφ); sin(Δφ) cos(Δφ)]
    d = inv(I(2) - R) * SA[x - x0, y - y0]
    return Point2f(x0 + d[1], y0 + d[2])
end

frp_pts = [frp_point(view(sol, :, i)) for i in axes(sol, 2)]

# Cross-check without any fit: the pole of the finite rotation of the moving link,
# taken from the solution itself, must coincide with the centroid point.
function mcs_deviation(x2, y2, φ2)
    x0, y0, φ0 = x2[1], y2[1], φ2[1]
    dev = 0.0
    for i in 2:length(φ2)-1
        dφ = φ2[i+1] - φ2[i-1]
        abs(dφ) < 1e-4 && continue

        ux = x2[i] - x0
        uy = y2[i] - y0
        dux = (x2[i+1] - x2[i-1]) / dφ
        duy = (y2[i+1] - y2[i-1]) / dφ
        icr = SA[x0 + ux - duy, y0 + uy + dux]      # МЦС из решения
        c = C1(φ2[i] - φ0)
        contact = SA[nominal[1] + c[1], nominal[2] + c[2]]
        dev = max(dev, hypot(icr[1] - contact[1], icr[2] - contact[2]))
    end
    return dev
end

dev = mcs_deviation(x2, y2, φ2)
println("МЦС из решения vs точка центроиды: max ", round(1e9 * dev; digits = 2), " nm")

##
# Two panels: the whole mechanism on the left and a close-up of the hinge on the
# right, where the travel of the МЦС along the centroids is clearly visible. The
# bounds are computed once from the solution, so the camera does not jitter while
# recording; the margin is wider than the default 5%, because the fixed-joint
# hatch sits above the link and `trajectory_limits` counts bodies only.
fig = Figure(size = (1200, 560))

ax = Axis(fig[1, 1], aspect = DataAspect())
ax.title = "Hinge with drift: whole mechanism"
p = systemplot!(ax, sys, sol)
limits!(ax, trajectory_limits(sys, sol; margin = 0.1)...)

axz = Axis(fig[1, 2], aspect = DataAspect())
axz.title = "hinge area: МЦС, ПКП and the centroids"
pz = systemplot!(axz, sys, sol)
limits!(axz, -0.05, 0.05, -link_len - 0.03, -link_len + 0.03)

mcs = Observable(mcs_pts[1])
frp = Observable(frp_pts[1])
moving_centroid = Observable(moving_centroid_pts(view(sol, :, 1)))
for axis in (ax, axz)
    lines!(axis, fixed_centroid_pts, color = :black, linewidth = 1.5,
           linestyle = :dash)
    lines!(axis, moving_centroid, color = :teal, linewidth = 1.5,
           linestyle = :dashdot)
    lines!(axis, frp_pts, color = :crimson, linewidth = 1.5)
    scatter!(axis, [nominal], color = :gray, marker = :x, markersize = 18)
    scatter!(axis, mcs, color = :orange, strokecolor = :black, strokewidth = 1,
             markersize = 15)
    scatter!(axis, frp, color = :transparent, strokecolor = :crimson, strokewidth = 2,
             markersize = 15)
end
# Подписи — на правой панели, оттуда же строится общая легенда.
lines!(axz, fixed_centroid_pts, color = :black, linewidth = 1.5,
       linestyle = :dash, label = "неподвижная центроида")
lines!(axz, moving_centroid, color = :teal, linewidth = 1.5,
       linestyle = :dashdot, label = "подвижная центроида")
lines!(axz, frp_pts, color = :crimson, linewidth = 1.5, label = "траектория ПКП")
scatter!(axz, mcs, color = :orange, strokecolor = :black, strokewidth = 1,
         markersize = 15, label = "МЦС")
scatter!(axz, frp, color = :transparent, strokecolor = :crimson, strokewidth = 2,
         markersize = 15, label = "ПКП")
scatter!(axz, [nominal], color = :gray, marker = :x, markersize = 18,
         label = "номинальный шарнир")
Legend(fig[2, 1:2], axz; orientation = :horizontal, framevisible = false,
       labelsize = 11, nbanks = 2)

fps = round(Int, 1 / (20 * step(time_span)))
record(fig, "out/hinge_with_drift.mp4", 1:20:length(time_span);
       framerate = fps) do i
    p.time = i
    pz.time = i
    mcs[] = mcs_pts[i]
    frp[] = frp_pts[i]
    moving_centroid[] = moving_centroid_pts(view(sol, :, i))
end

# Отдельный кадр максимального отклонения — для отчётов и документации.
imax = findfirst(t -> t >= time_span[end] / 4, time_span)
p.time = imax
pz.time = imax
mcs[] = mcs_pts[imax]
frp[] = frp_pts[imax]
moving_centroid[] = moving_centroid_pts(view(sol, :, imax))
frame_path = joinpath(@__DIR__, "..", "out", "hinge_with_drift.png")
save(frame_path, fig)
println("saved ", frame_path)

##
# Две сверки на графиках: МЦС, предписанный уравнениями центроид, против полюса,
# посчитанного по состояниям из симуляции (конечные разности), и ПКП
# аналитической позы против ПКП по состояниям.
t = collect(time_span)
θ_t = φ2 .- φ0                                   # угол дрейфа по кадрам
θ_cmd = [amplitude * sin(ω * ti) for ti in t]    # он же предписанный мотором

# МЦС по уравнениям центроид (аналитика) и по состояниям (полюс из разностей).
mcs_an = [C1(θ) for θ in θ_cmd]
mcs_sim = [SA[NaN, NaN] for _ in t]
for i in 2:length(t)-1
    dθ = θ_t[i+1] - θ_t[i-1]
    abs(dθ) < 1e-4 && continue

    dux = (x2[i+1] - x2[i-1]) / dθ
    duy = (y2[i+1] - y2[i-1]) / dθ
    ux = x2[i] - x0
    uy = y2[i] - y0
    mcs_sim[i] = SA[x0 + ux - duy - nominal[1], y0 + uy + dux - nominal[2]]
end

# ПКП: аналитическая поза из уравнений связи и поза из состояний.
R(φ) = SA[cos(φ) -sin(φ); sin(φ) cos(φ)]
θ1 = -π/2
pose2(θ) = SA[nominal[1], nominal[2]] + C1(θ) - R(θ1 + π + θ) * (h2 + c2(θ))
p20 = pose2(0.0)

function frp_from(p0, p, Δφ)
    abs(Δφ) < 1e-6 && return SA[NaN, NaN]       # вырождение: ПКП → МЦС
    return p0 + inv(I(2) - R(Δφ)) * (p - p0)
end

nominal_v = SA[nominal[1], nominal[2]]
frp_an = [frp_from(p20, pose2(θ), θ) - nominal_v for θ in θ_cmd]
frp_sim = [frp_from(SA[x0, y0], SA[x2[i], y2[i]], φ2[i] - φ0) - nominal_v
           for i in eachindex(t)]

dev_mcs = maximum(hypot(mcs_sim[i][1] - mcs_an[i][1], mcs_sim[i][2] - mcs_an[i][2])
                  for i in eachindex(t) if !isnan(mcs_sim[i][1]))
dev_frp = maximum(hypot(frp_sim[i][1] - frp_an[i][1], frp_sim[i][2] - frp_an[i][2])
                  for i in eachindex(t) if !isnan(frp_sim[i][1]))
println("МЦС: аналитика vs симуляция — max ", round(1e9 * dev_mcs; digits = 1), " нм")
println("ПКП: аналитика vs симуляция — max ", round(1e9 * dev_frp; digits = 1), " нм")

chk = Figure(size = (1200, 640))
ax11 = Axis(chk[1, 1], title = "МЦС: Δx от номинала", xlabel = "t, с", ylabel = "мм")
ax21 = Axis(chk[2, 1], title = "МЦС: Δy от номинала", xlabel = "t, с", ylabel = "мм")
ax12 = Axis(chk[1, 2], title = "ПКП: Δx от номинала", xlabel = "t, с", ylabel = "мм")
ax22 = Axis(chk[2, 2], title = "ПКП: Δy от номинала", xlabel = "t, с", ylabel = "мм")

for (axis, an, sim, k) in ((ax11, mcs_an, mcs_sim, 1), (ax21, mcs_an, mcs_sim, 2),
                           (ax12, frp_an, frp_sim, 1), (ax22, frp_an, frp_sim, 2))
    lines!(axis, t, 1e3 .* [v[k] for v in an], color = :black, linewidth = 4,
           alpha = 0.25, label = "аналитика")
    lines!(axis, t, 1e3 .* [v[k] for v in sim], color = :crimson, linewidth = 1.5,
           linestyle = :dash, label = "симуляция")
end
axislegend(ax11; position = :lt, framevisible = false, labelsize = 10)

check_path = joinpath(@__DIR__, "..", "out", "hinge_drift_check.png")
save(check_path, chk)
println("saved ", check_path)
