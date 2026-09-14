using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff
using StaticArrays

const g = 9.81

# ---------------------------------------------------------------------------
# Задаваемое движение ползуна вдоль направляющей:
#     s(t) = s0 + A*sin(Ω*t)
# Скорость и ускорение считаются аналитически и передаются в актуатор
# (иначе ∂f/∂t в якобиане будет не согласован с самим таргетом).
# ---------------------------------------------------------------------------
const s0 = 1.0
const A  = 0.5
const Ω  = 2.0

s_target(t) = s0 + A * sin(Ω * t)
v_target(t) = A * Ω * cos(Ω * t)
a_target(t) = -A * Ω^2 * sin(Ω * t)

# Аналитическая сила, которую должен развивать актуатор вдоль оси слайдера.
# Проецируем уравнение движения ползуна на единичную направляющую n⃗:
#     m * s̈ = λ + m * g⃗·n⃗ ,   g⃗ = (0, -g),  n⃗ = (cosα, sinα)
#   ⇒ λ_phys = m*s̈ + m*g*sinα
#
# ВНИМАНИЕ (результат проверки): текущая формулировка актуатора задаёт связь
# по положению (s = target), и при этом множитель Лагранжа, который выдаёт
# решатель CRos, равен НЕ физической силе, а
#     λ_cros = -m*s̈ + m*g*sinα
# (внешняя/гравитационная часть верна, а инерционная входит с обратным знаком).
# Положение s(t) и скорость v(t) при этом отслеживаются точно. Чтобы получить
# физическую силу λ, связь актуатора нужно задавать на уровне скоростей
# (index-1), а не положений.
λ_analytic(t, angle, mass) = mass * a_target(t) - mass * g * sin(angle)


# ---------------------------------------------------------------------------
# Сборка системы: bd1 зафиксировано, bd2 скользит по направляющей bd1
# ---------------------------------------------------------------------------
function build_system(angle_deg::Float64, time_span)
    angle = deg2rad(angle_deg)
    dir = SA[cos(angle), sin(angle)]

    bd1 = Body2D(1.0, 0.1)          # направляющая (зафиксирована)
    bd2 = Body2D(1.0, 0.1)          # ползун
    bd2.forces[2] = (state, t) -> -bd2.mass * g   # гравитация

    jnt1 = FixedJoint(bd1)
    setposition!(jnt1, SA[0.0, 0.0])
    setrotation!(jnt1, 0.0)

    jnt2 = SliderJoint(bd1, bd2)
    set_position_on_first_body!(jnt2, SA[0.0, 0.0])
    set_position_on_second_body!(jnt2, SA[0.0, 0.0])
    set_direction_on_first_body!(jnt2, dir)
    set_direction_on_second_body!(jnt2, dir)

    mot = PositionLinearActuator2D(jnt2, s_target(time_span[1]), v_target(time_span[1]))

    sys = MBSystem2D()
    add!(sys, bd1)
    add!(sys, bd2)
    add!(sys, jnt1)
    add!(sys, jnt2)
    add!(sys, mot)

    if (!assemble!(sys))
        error("Assembling failed!")
    end

    # Согласованное начальное состояние: центр ползуна = s(0)*n⃗, скорость = v(0)*n⃗
    initial = zeros(number_of_dofs(sys))
    bd2_p = get_body_position_dofs(sys, bd2)
    bd2_v = get_body_velocity_dofs(sys, bd2)
    initial[bd2_p[1]] = s_target(time_span[1]) * dir[1]
    initial[bd2_p[2]] = s_target(time_span[1]) * dir[2]
    initial[bd2_v[1]] = v_target(time_span[1]) * dir[1]
    initial[bd2_v[2]] = v_target(time_span[1]) * dir[2]

    sys.prestep = (state) -> begin
        t = state[end]
        settarget!(mot, s_target(t), v_target(t))
    end

    return sys, initial, mot, bd2, angle
end

# ---------------------------------------------------------------------------
# Прогоняем горизонтальный и наклонный слайдер и сравниваем λ с аналитикой
# ---------------------------------------------------------------------------
time_span = 0:0.0002:6.0
cases = [0.0, 30.0]

results = Any[]

for angle_deg in cases
    sys, initial, mot, bd2, angle = build_system(angle_deg, time_span)
    sol = simulate(sys, initial, time_span)

    λ_sim = sol[get_lms(sys, mot)[1], :]
    λ_exact = λ_analytic.(time_span, angle, 1.0)
    λ_cros_exact = λ_cros.(time_span, angle, 1.0)

    # Скорость ползуна, спроецированная на направляющую
    bd2_v = get_body_velocity_dofs(sys, bd2)
    v_sim = @. sol[bd2_v[1], :] * cos(angle) + sol[bd2_v[2], :] * sin(angle)
    v_exact = v_target.(time_span)

    # Координата ползуна вдоль направляющей
    bd2_p = get_body_position_dofs(sys, bd2)
    s_sim = @. sol[bd2_p[1], :] * cos(angle) + sol[bd2_p[2], :] * sin(angle)
    s_exact = s_target.(time_span)

    # Множители в начальном состоянии нулевые, поэтому на первом шаге есть
    # короткий переходный процесс — силовые ошибки считаем при t > 0.05.
    tsel = time_span .> 0.05

    println("angle = $(angle_deg)°")
    println("  max |s - s_target|        = $(round(maximum(abs.(s_sim .- s_exact)), digits = 8))")
    println("  max |v·n - v_target|      = $(round(maximum(abs.(v_sim .- v_exact)), digits = 8))")
    println("  max |λ - (m·s̈ + m·g·sinα)| = $(round(maximum(abs.(λ_sim[tsel] .- λ_exact[tsel])), digits = 8))   <- физическая сила")


    push!(results, (angle_deg = angle_deg, sys = sys, sol = sol, λ_num = λ_sim,
                    λ_exact = λ_exact))
end

# ==== PLOTTING ====
##
fig = Figure(size = (1100, 500))
for (i, res) in enumerate(results)
    ax = Axis(fig[i, 1], title = "Слайдер $(res.angle_deg)°",
              xlabel = "t, с", ylabel = "λ, Н")
    lines!(ax, time_span, res.λ_num - res.λ_exact, label = "error")
    # lines!(ax, time_span, res.λ_num, label = "λ из симуляции")
    # lines!(ax, time_span, res.λ_exact, linestyle = :dash, label = "m·s̈ + m·g·sinα (физическая)")
    axislegend(ax)
end
fig
##
# Визуализация наклонного случая
last = results[end]
animate(last.sys, last.sol, time_span, "out/linear_actuator.mp4";
        framerate = 30, limits = (-0.5, 2.0, -0.5, 1.5))

