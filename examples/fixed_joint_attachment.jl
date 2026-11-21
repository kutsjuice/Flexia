using Pkg; Pkg.activate(@__DIR__)
using Flexia
using CairoMakie
using StaticArrays

# ===========================================================================
# FixedVisual + Attachment: где и как лежит штриховка заделки
# ===========================================================================
#
# Самая маленькая модель: одно тело и одна заделка. Динамики нет — состояние
# одно, кадр один, поэтому видно только геометрию.
#
# Панели отличаются ровно одним: `attachment` у `FixedVisual`. По умолчанию
# полоса стоит под звеном; `offset` сдвигает её в системе тела, `angle`
# доворачивает её относительно тела. Длина полосы — длина звена, высота —
# четверть длины, поэтому у звена длиной 2 полоса имеет размер 2 × 0.5.
#
# Запуск (из корня пакета, там же появится out/fixed_joint_attachment.png):
#
#     julia --project=examples examples/fixed_joint_attachment.jl

"""
    fixed_joint_case(att::Attachment) -> (sys, sol)

Одно тело `2 × 0.3` в начале координат и заделка на нём с геометрией `att`.
Готовое «решение» из одного кадра — его достаточно для `systemplot!`.
"""
function fixed_joint_case(att::Attachment)
    bd = Body2D(1.0, 1.0; visual = BarVisual(length = 2.0, width = 0.3))

    sys = MBSystem2D()
    add!(sys, bd)

    fj = FixedJoint(bd, visual = FixedVisual(attachment = att))
    add!(sys, fj)

    setposition!(fj, SA[0.0, 0.0])
    setrotation!(fj, 0.0)

    assemble!(sys) || error("assemble! failed")

    state = zeros(number_of_dofs(sys))
    location = get_body_position_dofs(sys, bd)
    state[location[1]], state[location[2]], state[location[3]] = 0.0, 0.0, 0.0

    return sys, reshape(state, :, 1)
end

cases = [
    ("Attachment()", Attachment()),
    ("offset = (0, -0.9)", Attachment(offset = SA[0.0, -0.9])),
    ("angle = -π/2", Attachment(angle = -π / 2)),
    ("offset + angle = π/4", Attachment(offset = SA[1.0, -0.9], angle = π / 4)),
]

fig = Figure(size = (1000, 820))
Label(fig[0, 1:2], "Меняется только FixedVisual.attachment (полоса заделки)",
      fontsize = 16, halign = :center)

for (i, (title, att)) in enumerate(cases)
    sys, sol = fixed_joint_case(att)

    row, col = (i - 1) ÷ 2 + 1, (i - 1) % 2 + 1
    ax = Axis(fig[row, col], aspect = DataAspect(), title = title,
              xticks = -3:1:3, yticks = -3:1:3)
    systemplot!(ax, sys, sol)
    # Границы одинаковые у всех панелей, иначе масштаб мешает сравнивать.
    limits!(ax, -3.0, 3.0, -3.0, 3.0)
end

mkpath("out")
save("out/fixed_joint_attachment.png", fig)
println("saved out/fixed_joint_attachment.png")
