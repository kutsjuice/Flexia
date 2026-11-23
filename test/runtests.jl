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

        # --- Attachment полосы: смещение и доворот в системе тела -----------
        # Отдельная система «тело + заделка», чтобы не менять общий фикстур:
        # его заделка участвует в тесте рендера ниже.
        function hatch_with(att; θ = 0.0, offset = SA[0.0, 0.0])
            bd = Body2D(1.0, 1.0; length = l)
            sys2 = MBSystem2D()
            add!(sys2, bd)
            fj2 = FixedJoint(bd, visual = FixedVisual(attachment = att))
            add!(sys2, fj2)
            setposition!(fj2, SA[offset[1], offset[2]])
            setrotation!(fj2, θ)
            assemble!(sys2) || error("assemble! failed")
            st = zeros(number_of_dofs(sys2))
            d = get_body_position_dofs(sys2, bd)
            st[d[1]], st[d[2]], st[d[3]] = offset[1], offset[2], θ
            return Flexia.fixed_joint_hatch_points(sys2, fj2, st)
        end

        # поворот тела переносит полосу вместе с ним: при θ = π/2 t-ось
        # «вниз, под звеном» уходит в +x
        rot = hatch_with(Attachment(); θ = π / 2)
        @test all(-1e-6 <= p[1] <= h + 1e-6 for p in rot)
        @test all(-l / 2 - 1e-6 <= p[2] <= l / 2 + 1e-6 for p in rot)

        # смещение сдвигает полосу целиком
        off = hatch_with(Attachment(offset = SA[1.0, 2.0]))
        @test all(1.0 - l / 2 - 1e-6 <= p[1] <= 1.0 + l / 2 + 1e-6 for p in off)
        @test all(2.0 - h - 1e-6 <= p[2] <= 2.0 + 1e-6 for p in off)

        # доворот на -π/2 разворачивает полосу относительно звена
        ang = hatch_with(Attachment(angle = -π / 2))
        @test all(-h - 1e-6 <= p[1] <= 1e-6 for p in ang)
        @test all(-l / 2 - 1e-6 <= p[2] <= l / 2 + 1e-6 for p in ang)
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
    # Шарнир с дрейфом: точка контакта уходит от номинальной вдоль
    # аналитических центроид, относительный угол отсчитывается от rot_offset.
    @testset "шарнир с дрейфом (МЦС вдоль центроид)" begin
        a, b = 0.4, 0.15                 # точка крепления на теле 2
        k, c = 0.02, 0.005               # наклоны центроид, м/рад
        drift1 = θ -> SA[k * θ, 0.0]     # неподвижная центроида на теле 1
        drift2 = θ -> SA[0.0, c * θ]     # подвижная центроида на теле 2

        function hinge_pair(; d1 = nothing, d2 = nothing, offset = 0.0)
            bd1 = Body2D(1.0, 1.0; length = 2.0)
            bd2 = Body2D(1.0, 1.0; length = 1.0)
            s = MBSystem2D()
            add!(s, bd1)
            add!(s, bd2)
            hj = HingeJoint(bd1, bd2; drift1 = d1, drift2 = d2, rot_offset = offset)
            set_position_on_first_body!(hj, SA[0.0, 0.0])
            set_position_on_second_body!(hj, SA[a, b])
            add!(s, hj)
            assemble!(s)
            return s, hj, bd1, bd2
        end

        function hinge_rhs(s, hj, st)
            r = zeros(number_of_dofs(s))
            Flexia.add_to_rhs!(r, st, s, hj)
            return r
        end

        # Прежняя формула (только номинальные точки) — эталон для проверки
        # обратной совместимости.
        function legacy_hinge_rhs!(r, st, s, hj)
            bd1 = hj.body1
            bd2 = hj.body2
            q1 = get_body_position_dofs(s, bd1)
            u1 = get_body_velocity_dofs(s, bd1)
            q2 = get_body_position_dofs(s, bd2)
            u2 = get_body_velocity_dofs(s, bd2)

            xi, yi, θi = st[q1]
            xj, yj, θj = st[q2]

            lms = get_lms(s, hj)
            λ1 = st[lms[1]]
            λ2 = st[lms[2]]

            xci = hj.body1_hinge_point[1]
            yci = hj.body1_hinge_point[2]
            xcj = hj.body2_hinge_point[1]
            ycj = hj.body2_hinge_point[2]

            r[u1[1]] += λ1
            r[u1[2]] += λ2
            r[u1[3]] += (-xci * sin(θi) - yci * cos(θi)) * λ1 +
                        (xci * cos(θi) - yci * sin(θi)) * λ2

            r[u2[1]] += -λ1
            r[u2[2]] += -λ2
            r[u2[3]] += (xcj * sin(θj) + ycj * cos(θj)) * λ1 +
                        (-xcj * cos(θj) + ycj * sin(θj)) * λ2

            r[lms[1]] = (xi + xci * cos(θi) - yci * sin(θi)) -
                        (xj + xcj * cos(θj) - ycj * sin(θj))
            r[lms[2]] = (yi + xci * sin(θi) + yci * cos(θi)) -
                        (yj + xcj * sin(θj) + ycj * cos(θj))
            return nothing
        end

        # --- обратная совместимость: обычный шарнир ------------------------
        sys0, hj0, bd10, bd20 = hinge_pair()
        @test hj0 isa HingeJoint{Nothing, Nothing}     # параметр типа не мешает
        @test number_of_dofs(hj0) == 2
        @test hj0.rot_offset == 0.0
        @test visual(hj0) isa HingeVisual
        @test Flexia._eval_drift(hj0.drift1, 1.7) == SA[0.0, 0.0]
        @test Flexia._eval_drift(hj0.drift2, -0.3) == SA[0.0, 0.0]

        q1 = get_body_position_dofs(sys0, bd10)
        q2 = get_body_position_dofs(sys0, bd20)
        lms0 = get_lms(sys0, hj0)

        st0 = zeros(number_of_dofs(sys0))
        st0[q1] .= SA[0.3, -0.2, 0.4]
        st0[q2] .= SA[1.1, 0.25, -0.15]
        st0[lms0[1]] = 1.5
        st0[lms0[2]] = -0.7
        st0[end] = 1.0

        new0 = zeros(number_of_dofs(sys0))
        old0 = zeros(number_of_dofs(sys0))
        Flexia.add_to_rhs!(new0, st0, sys0, hj0)
        legacy_hinge_rhs!(old0, st0, sys0, hj0)
        @test isapprox(new0, old0; rtol = 1e-15, atol = 0)

        # дрейф не добавляет аллокаций обычному шарниру
        Flexia.add_to_rhs!(new0, st0, sys0, hj0)       # прогрев
        @test (@allocated Flexia.add_to_rhs!(new0, st0, sys0, hj0)) == 0

        # --- rot_offset: drift получает θ = θj - θi - offset ---------------
        seen = Ref(NaN)
        spy = θ -> begin
            seen[] = θ
            return SA[0.0, 0.0]
        end
        sys_r, hj_r, bd1r, bd2r = hinge_pair(d1 = spy, offset = 0.25)
        st_r = zeros(number_of_dofs(sys_r))
        st_r[get_body_position_dofs(sys_r, bd1r)] .= SA[0.0, 0.0, 0.0]
        st_r[get_body_position_dofs(sys_r, bd2r)] .= SA[1.0, 0.0, 0.25]
        st_r[end] = 1.0

        hinge_rhs(sys_r, hj_r, st_r)
        @test seen[] === 0.0        # θj - θi == rot_offset → ровно ноль
        @test hj_r.rot_offset == 0.25

        set_rot_offset!(hj_r, 0.0)
        @test hj_r.rot_offset == 0.0
        hinge_rhs(sys_r, hj_r, st_r)
        @test seen[] === 0.25

        # --- аналитическая центроида: поворот вокруг текущей точки МЦС -----
        sys_d, hj_d, bd1d, bd2d = hinge_pair(d1 = drift1, d2 = drift2)
        @test hj_d isa HingeJoint{typeof(drift1), typeof(drift2)}

        θrel = 0.25
        xi, yi, θi = 0.3, -0.2, 0.1
        θj = θi + θrel

        rot(φ) = SA[cos(φ) -sin(φ); sin(φ) cos(φ)]
        w = SA[xi, yi] + rot(θi) * drift1(θrel)              # точка МЦС в мире
        p2 = w - rot(θj) * (SA[a, b] + drift2(θrel))         # поза тела 2

        q1d = get_body_position_dofs(sys_d, bd1d)
        q2d = get_body_position_dofs(sys_d, bd2d)
        lmsd = get_lms(sys_d, hj_d)

        st_d = zeros(number_of_dofs(sys_d))
        st_d[q1d] .= SA[xi, yi, θi]
        st_d[q2d] .= SA[p2[1], p2[2], θj]
        st_d[lmsd[1]] = 1.5
        st_d[lmsd[2]] = -0.7
        st_d[end] = 1.0

        rhs_d = hinge_rhs(sys_d, hj_d, st_d)
        @test abs(rhs_d[lmsd[1]]) < 1e-15      # невязка шарнира = 0
        @test abs(rhs_d[lmsd[2]]) < 1e-15

        # проверка не «пустая»: сдвиг тела 2 ломает связь
        perturbed = copy(st_d)
        perturbed[q2d[2]] += 1e-3
        @test hinge_rhs(sys_d, hj_d, perturbed)[lmsd[2]] ≈ -1e-3

        # без дрейфа та же поза шарниру уже не удовлетворяет
        sys_nd, hj_nd = hinge_pair()
        @test abs(hinge_rhs(sys_nd, hj_nd, st_d)[lmsd[1]]) > 1e-3

        # маркер шарнира стоит на текущей точке МЦС, а не на номинальной
        @test Flexia.get_hinge_point(sys_d, hj_d, st_d) ≈ Point2f(w) atol = 1e-6

        # --- якобиан: Dual доходит до drift-функций без ошибок типов -------
        n = number_of_dofs(sys_d)
        J = sys_d.jacobian(st_d)
        @test size(J) == (n, n)
        @test all(isfinite, J)
        # дрейф виден в строке невязки: производная по углу тела 2 не нулевая
        @test abs(J[lmsd[1], q2d[3]]) > 0
        @test !isapprox(J[lmsd, :], sys_nd.jacobian(st_d)[lmsd, :])
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
