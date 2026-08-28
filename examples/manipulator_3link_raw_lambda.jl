using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff
using StaticArrays

const g = 9.81

# ---------------------------------------------------------------------------
# 1. Bodies: fixed ground + 3 rigid links (base -> shoulder -> elbow -> wrist)
# ---------------------------------------------------------------------------
L1, L2, L3 = 0.40, 0.35, 0.25                 # link lengths, m
m1, m2, m3 = 1.00, 0.80, 0.50                 # link masses, kg
I1, I2, I3 = m1 * L1^2 / 12, m2 * L2^2 / 12, m3 * L3^2 / 12   # slender-rod inertia

ground = Body2D(1e6, 1e6; length=0.05)        # "infinite" mass/inertia -> effectively fixed
link1  = Body2D(m1, I1; length=L1)
link2  = Body2D(m2, I2; length=L2)
link3  = Body2D(m3, I3; length=L3)

link1.forces[2] = (x, t) -> -link1.mass * g
link2.forces[2] = (x, t) -> -link2.mass * g
link3.forces[2] = (x, t) -> -link3.mass * g

# ---------------------------------------------------------------------------
# 2. Joints: ground is pinned to the world, links are chained by hinges
# ---------------------------------------------------------------------------
ground_joint = FixedJoint(ground)
setposition!(ground_joint, SA[0.0, 0.0])
setrotation!(ground_joint, 0.0)

hinge1 = HingeJoint(ground, link1)             # base joint
set_position_on_first_body!(hinge1, SA[0.0, 0.0])
set_position_on_second_body!(hinge1, SA[-L1/2, 0.0])

hinge2 = HingeJoint(link1, link2)              # shoulder/elbow joint
set_position_on_first_body!(hinge2, SA[L1/2, 0.0])
set_position_on_second_body!(hinge2, SA[-L2/2, 0.0])

hinge3 = HingeJoint(link2, link3)              # wrist joint
set_position_on_first_body!(hinge3, SA[L2/2, 0.0])
set_position_on_second_body!(hinge3, SA[-L3/2, 0.0])

motor1 = PositionMotor2D(hinge1, 0.0)
motor2 = PositionMotor2D(hinge2, 0.0)
motor3 = PositionMotor2D(hinge3, 0.0)

sys = MBSystem2D()
add!(sys, ground); add!(sys, link1); add!(sys, link2); add!(sys, link3)
add!(sys, ground_joint); add!(sys, hinge1); add!(sys, hinge2); add!(sys, hinge3)
add!(sys, motor1); add!(sys, motor2); add!(sys, motor3)

if !assemble!(sys)
    println("Assembling failed!")
end

# ---------------------------------------------------------------------------
# 3. Initial configuration: arm fully extended along +x (q1=q2=q3=0)
# ---------------------------------------------------------------------------
q0 = (0.0, 0.0, 0.0)          # initial joint angles, rad
qf = (π/2, -π/3, -π/4)        # target joint angles, rad ("reach up and fold")

initial_state = zeros(number_of_dofs(sys))

θ1_0 = q0[1]
θ2_0 = θ1_0 + q0[2]
θ3_0 = θ2_0 + q0[3]

x1, y1 = (L1/2) * cos(θ1_0), (L1/2) * sin(θ1_0)
ex1, ey1 = L1 * cos(θ1_0), L1 * sin(θ1_0)

x2, y2 = ex1 + (L2/2) * cos(θ2_0), ey1 + (L2/2) * sin(θ2_0)
ex2, ey2 = ex1 + L2 * cos(θ2_0), ey1 + L2 * sin(θ2_0)

x3, y3 = ex2 + (L3/2) * cos(θ3_0), ey2 + (L3/2) * sin(θ3_0)

set_initial_position!(initial_state, sys, link1, SA_F64[x1, y1, θ1_0])
set_initial_position!(initial_state, sys, link2, SA_F64[x2, y2, θ2_0])
set_initial_position!(initial_state, sys, link3, SA_F64[x3, y3, θ3_0])

# ---------------------------------------------------------------------------
# 4. Trajectory: smooth (quintic, zero vel/accel at the ends) point-to-point
#    move in joint space, executed over T_move seconds, then held.
# ---------------------------------------------------------------------------
T_move = 2.0
quintic(τ) = (τ = clamp(τ, 0.0, 1.0); 10τ^3 - 15τ^4 + 6τ^5)

sys.prestep = (state) -> begin
    t = state[end]
    s = quintic(t / T_move)
    settarget!(motor1, q0[1] + (qf[1] - q0[1]) * s)
    settarget!(motor2, q0[2] + (qf[2] - q0[2]) * s)
    settarget!(motor3, q0[3] + (qf[3] - q0[3]) * s)
end

# ---------------------------------------------------------------------------
# 5. Simulate (2 s of motion + 2 s holding the final pose)
# ---------------------------------------------------------------------------
time_span = 0:0.001:4.0
sol = simulate(sys, initial_state, time_span)

animate(sys, sol, time_span, "out/manipulator_3link_raw.mp4";
    framerate=floor(Int64, 1.0 / step(time_span)), limits=(-1, 1, -1, 1))

# ---------------------------------------------------------------------------
# 6. Plot the λ rows straight out of `sol`, with no post-processing. This is
#    the multiplier `cros!` itself carries in the state - cheap to get, but
#    per `honest_constraint_forces`'s docstring it is only an index-1
#    (velocity-consistent) reconstruction of a constraint enforced at the
#    position level, so it is known to drift away from the true constraint
#    force as the time step shrinks (no Baumgarte/GGL stabilization here).
#    Compare against `manipulator_3link.jl`, which reconstructs the same
#    joint torques via `honest_constraint_forces` instead.
# ---------------------------------------------------------------------------
τ1 = sol[get_lms(sys, motor1)[1], :]
τ2 = sol[get_lms(sys, motor2)[1], :]
τ3 = sol[get_lms(sys, motor3)[1], :]

fig_torques = Figure(size=(900, 500))
ax_t = Axis(fig_torques[1, 1], xlabel="t, s", ylabel="τ, N·m",
    title="Момент в приводах (сырые λ из sol)")
lines!(ax_t, time_span, τ1, label="Сустав 1 (основание)")
lines!(ax_t, time_span, τ2, label="Сустав 2 (плечо)")
lines!(ax_t, time_span, τ3, label="Сустав 3 (кисть)")
vlines!(ax_t, [T_move], color=:gray, linestyle=:dash, label="конец движения")
axislegend(ax_t, position=:rb)
save("out/manipulator_3link_raw_torques.png", fig_torques)

# Full breakdown: every λ component of every connector, straight from `sol`.
ncols = 3
nrows = cld(length(sys.connectors), ncols)
fig_all = Figure(size=(1300, 320 * nrows))
for (i, conn) in enumerate(sys.connectors)
    r, c = fldmod1(i, ncols)
    ax = Axis(fig_all[r, c], xlabel="t, s", ylabel="λ",
        title=string(nameof(typeof(conn))) * " #$i")
    for (j, idx) in enumerate(get_lms(sys, conn))
        lines!(ax, time_span, sol[idx, :], label="λ$j")
    end
    axislegend(ax, position=:rb, labelsize=10)
end
save("out/manipulator_3link_raw_all_lambdas.png", fig_all)

fig_torques
