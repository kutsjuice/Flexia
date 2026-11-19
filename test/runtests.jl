using Flexia
using Test
using StaticArrays
using LinearAlgebra
using CairoMakie          # тестовая зависимость: рендер без дисплея

# ---------------------------------------------------------------------------
# Пользовательская форма тела: проверяем, что новая геометрия добавляется
# БЕЗ правок математических файлов и существующих рецептов.
# ---------------------------------------------------------------------------
struct HexVisual <: AbstractBodyVisual
    radius::Float64
end

function Flexia.body_contour(x, y, θ, v::HexVisual)
    att = Attachment()
    return Point2f[
        Flexia.to_global(x, y, θ, att,
                         Point2f(v.radius * cos(a), v.radius * sin(a)))
        for a in range(0, 2π; length = 7)
    ]
end

function Flexia._draw_body!(p::Flexia.BodyPlot, v::HexVisual)
    pts = Makie.lift(p.time) do t
        Flexia.body_contour(p.system[], p.body[], view(p.solution[], :, t))
    end
    Makie.lines!(p, pts; color = :black)
    return p
end

# Сущность вообще без визуального компаньона: не рисуется и не требует
# никаких методов отрисовки.
struct UnvisualizedEntity <: AbstractConnector2D end

# ---------------------------------------------------------------------------
# Общая модель для тестов геометрии
# ---------------------------------------------------------------------------
function geometry_fixture()
    bd1 = Body2D(1.0, 1.0; length = 2.0)
    bd2 = Body2D(1.0, 1.0; visual = RectVisual(width = 0.4, height = 0.3))
    bd3 = Body2D(1.0, 1.0; visual = CircleVisual(radius = 0.25))

    sys = MBSystem2D()
    for bd in (bd1, bd2, bd3)
        add!(sys, bd)
    end

    fj = FixedJoint(bd1)
    hj = HingeJoint(bd1, bd2)
    sj = SliderJoint(bd1, bd3)
    for cn in (fj, hj, sj)
        add!(sys, cn)
    end
    ls = LinearSpring(sj; vis_r = 0.2, vis_N = 6)
    add!(sys, ls)
    assemble!(sys) || error("assemble! failed")

    state = zeros(number_of_dofs(sys))
    state[get_body_position_dofs(sys, bd1)] .= SA[0.0, 0.0, 0.0]
    state[get_body_position_dofs(sys, bd2)] .= SA[1.0, 0.5, 0.3]
    state[get_body_position_dofs(sys, bd3)] .= SA[0.5, 0.3, 0.4]
    state[end] = 1.0

    return sys, state, (; bd1, bd2, bd3, fj, hj, sj, ls)
end

@testset "Flexia.jl" begin
    @test test_func() == 1

    sys, state, e = geometry_fixture()

    # -----------------------------------------------------------------------
    @testset "визуальные компаньоны" begin
        @test visual(e.bd1) isa BarVisual
        @test visual(e.bd1).length == 2.0                 # kwarg length= проброшен
        @test visual(e.bd2) isa RectVisual
        @test visual(e.bd3) isa CircleVisual
        @test visual(e.fj) isa FixedVisual
        @test visual(e.hj) isa HingeVisual
        @test visual(e.sj) isa SliderVisual
        @test visual(e.sj).axis_length == 2.0             # значение по умолчанию
        @test visual(e.ls) isa LinearSpringVisual
        @test visual(e.ls).coil_radius == 0.2             # kwarg vis_r= проброшен
        @test visual(e.ls).n_coils == 6                   # kwarg vis_N= проброшен

        # цвет по умолчанию — палитра Makie
        @test visual(e.bd1).style.color === nothing
        # у заделки цвет задан явно (исторический первый цвет Wong)
        @test visual(e.fj).style.color == Makie.wong_colors()[1]

        # сущность без компаньона: NoVisual, не рисуется, методов не требует
        u = UnvisualizedEntity()
        @test visual(u) isa NoVisual
        @test !Flexia.is_drawn(u)
        @test !Flexia.is_drawn_visual(NoVisual())

        # скрыть конкретный экземпляр, не меняя тип
        @test Flexia.is_drawn(e.bd1)
        set_visualization!(e.bd1, BarVisual(length = 2.0, style = VisualStyle(visible = false)))
        @test !Flexia.is_drawn(e.bd1)
        set_visualization!(e.bd1, BarVisual(length = 2.0))
        @test Flexia.is_drawn(e.bd1)

        # пользовательский визуал без поля style не ломает ни видимость, ни zorder
        hex = HexVisual(0.2)
        @test Flexia.is_drawn_visual(hex)
        @test Flexia.style_of(hex) isa VisualStyle
    end

    # -----------------------------------------------------------------------
    @testset "контуры тел" begin
        bar = Flexia.body_contour(sys, e.bd1, state)
        @test length(bar) == 2
        @test bar[1] ≈ Point2f(-1.0, 0.0)
        @test bar[2] ≈ Point2f(1.0, 0.0)

        rect = Flexia.body_contour(sys, e.bd2, state)
        @test length(rect) == 5
        @test rect[1] ≈ rect[end]                          # контур замкнут
        @test norm(rect[2] - rect[1]) ≈ 0.4 atol = 1e-6    # размеры сохраняются
        @test norm(rect[3] - rect[2]) ≈ 0.3 atol = 1e-6    # при повороте

        circ = Flexia.body_contour(sys, e.bd3, state)
        @test length(circ) == 33
        @test circ[1] ≈ circ[end]
        @test all(isapprox(norm(p - Point2f(0.5, 0.3)), 0.25; atol = 1e-6) for p in circ)

        # пользовательская форма работает через тот же интерфейс
        hex = Flexia.body_contour(0.0, 0.0, 0.0, HexVisual(0.2))
        @test length(hex) == 7
        @test all(isapprox(norm(p), 0.2; atol = 1e-6) for p in hex)

        # NoVisual — геометрии нет
        @test isempty(Flexia.body_contour(0.0, 0.0, 0.0, NoVisual()))
    end

    # -----------------------------------------------------------------------
    @testset "штриховка заделки" begin
        pts = Flexia.fixed_joint_hatch_points(sys, e.fj, state)
        @test !isempty(pts)
        @test iseven(length(pts))                          # плоский список пар точек

        l = visual(e.bd1).length
        h = l * visual(e.fj).hatch_height_ratio
        # при θ = 0 локальные оси совпадают с глобальными: s = x, t = -y
        for p in pts
            @test -l / 2 - 1e-9 <= p[1] <= l / 2 + 1e-9
            @test -1e-9 <= -p[2] <= h + 1e-9
        end
    end

    # -----------------------------------------------------------------------
    @testset "ползун" begin
        axis_pts = Flexia.slider_axis_points(sys, e.sj, state)
        @test length(axis_pts) == 2
        @test norm(axis_pts[2] - axis_pts[1]) ≈ visual(e.sj).axis_length

        prism = Flexia.slider_prism_points(sys, e.sj, state)
        @test length(prism) == 5
        @test prism[1] ≈ prism[end]
        @test norm(prism[2] - prism[1]) ≈ visual(e.sj).slider_length atol = 1e-6
        @test norm(prism[3] - prism[2]) ≈ visual(e.sj).slider_width atol = 1e-6
    end

    # -----------------------------------------------------------------------
    @testset "линейная пружина" begin
        pts = Flexia.linear_spring_points(sys, e.ls, state)
        @test length(pts) == visual(e.ls).n_coils
        # первый конец — точно в точке крепления на теле 1
        d1 = get_body_position_dofs(sys, e.sj.body1)
        p1 = Point2f(state[d1[1]] + e.sj.body1_position[1],
                     state[d1[2]] + e.sj.body1_position[2])
        @test pts[1] ≈ p1
    end

    # -----------------------------------------------------------------------
    @testset "автолимиты" begin
        sol = repeat(state, 1, 5)
        lim = trajectory_limits(sys, sol)
        @test length(lim) == 4
        xmin, xmax, ymin, ymax = lim
        @test xmin < xmax
        @test ymin < ymax
        for t in 1:5, bd in bodies(sys)
            Flexia.is_drawn(bd) || continue
            for p in Flexia.body_contour(sys, bd, view(sol, :, t))
                @test xmin <= p[1] <= xmax
                @test ymin <= p[2] <= ymax
            end
        end

        # вырожденный случай: всё в одной точке — кадр не должен схлопнуться
        sys0 = MBSystem2D()
        bd0 = Body2D(1.0, 1.0; length = 0.0)
        add!(sys0, bd0)
        assemble!(sys0)
        lim0 = trajectory_limits(sys0, zeros(number_of_dofs(sys0), 1))
        @test lim0[1] < lim0[2]
        @test lim0[3] < lim0[4]
    end

    # -----------------------------------------------------------------------
    @testset "звено со скруглёнными концами (Link2DVisual)" begin
        l, w = 1.0, 0.2
        r = w / 2
        v = Link2DVisual(length = l, width = w)
        pts = Flexia.body_contour(0.0, 0.0, 0.0, v)

        @test pts[1] ≈ pts[end]                          # контур замкнут
        # 2 (низ) + 16 (правая дуга) + 1 (верх) + 16 (левая дуга)
        @test length(pts) == 35

        # габарит = прямоугольник плюс два скругления по width/2
        xs = first.(pts)
        ys = last.(pts)
        @test maximum(xs) ≈ l / 2 + r atol = 1e-6
        @test minimum(xs) ≈ -(l / 2 + r) atol = 1e-6
        @test maximum(ys) ≈ r atol = 1e-6
        @test minimum(ys) ≈ -r atol = 1e-6

        # все точки лежат в объединении прямоугольника и двух кружков.
        # Допуск 1e-6, а не 1e-9: Point2f — это Float32, и 0.1f0 не совпадает
        # с 0.1 в пределах 1e-9.
        tol = 1e-6
        inside(p) = (abs(p[1]) <= l / 2 + tol && abs(p[2]) <= r + tol) ||
                    hypot(p[1] - l / 2, p[2]) <= r + tol ||
                    hypot(p[1] + l / 2, p[2]) <= r + tol
        @test all(inside(p) for p in pts)

        # точки за прямоугольником лежат ровно на окружностях радиуса width/2
        right_arc = [p for p in pts if p[1] > l / 2 + tol]
        left_arc = [p for p in pts if p[1] < -l / 2 - tol]
        # по 16 точек на дугу, но крайние (в углах прямоугольника) сюда не попадают
        @test length(right_arc) == 15
        @test length(left_arc) == 15
        @test all(isapprox(hypot(p[1] - l / 2, p[2]), r; atol = tol) for p in right_arc)
        @test all(isapprox(hypot(p[1] + l / 2, p[2]), r; atol = tol) for p in left_arc)
        # в углах прямоугольника дуга касается граней: y = ±r
        @test any(isapprox(p[2], r; atol = tol) && isapprox(p[1], l / 2; atol = tol) for p in pts)
        @test any(isapprox(p[2], -r; atol = tol) && isapprox(p[1], -l / 2; atol = tol) for p in pts)

        # поворот тела переносит габарит в другую ось
        rot = Flexia.body_contour(0.0, 0.0, π / 2, v)
        @test maximum(last.(rot)) ≈ l / 2 + r atol = 1e-6
        @test maximum(abs.(first.(rot))) ≈ r atol = 1e-6

        # смещение из Attachment учитывается
        off = Flexia.body_contour(0.0, 0.0, 0.0,
                                  Link2DVisual(length = l, width = w,
                                         attachment = Attachment(offset = SA[1.0, 2.0])))
        @test maximum(first.(off)) ≈ 1.0 + l / 2 + r atol = 1e-6
        @test maximum(last.(off)) ≈ 2.0 + r atol = 1e-6
    end

    # -----------------------------------------------------------------------
    @testset "рендер без дисплея (CairoMakie)" begin
        sol = repeat(state, 1, 3)
        fig, ax, p = systemfigure(sys, sol)
        @test p.time[] == 1

        p.time = 2
        @test p.time[] == 2

        mktempdir() do dir
            path = joinpath(dir, "frame.png")
            save(path, fig)
            @test isfile(path)
            @test filesize(path) > 0
        end

        # новая форма тела рисуется через композитный рецепт без правок ядра
        sys2 = MBSystem2D()
        hexbd = Body2D(1.0, 1.0; visual = HexVisual(0.3))
        add!(sys2, hexbd)
        assemble!(sys2)
        sol2 = zeros(number_of_dofs(sys2), 2)
        fig2, _, p2 = systemfigure(sys2, sol2)
        @test length(p2.plots) == 1
        mktempdir() do dir
            path = joinpath(dir, "hex.png")
            save(path, fig2)
            @test isfile(path) && filesize(path) > 0
        end
    end
end
