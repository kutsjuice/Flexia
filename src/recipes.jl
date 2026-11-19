# ===========================================================================
# Рецепты Makie
# ===========================================================================
#
# Каждая сущность рисуется своим рецептом; рецепт выбирается по типу визуального
# компаньона. Композитный `SystemPlot` собирает всё в одну ось.
#
# Зависимости от бэкенда нет: файл использует только `Makie`. Бэкенд
# (GLMakie / CairoMakie / WGLMakie) выбирает скрипт пользователя через `using`.
#
# ПРО ЦВЕТ. `companion.style.color === nothing` означает «взять цвет из палитры
# Makie». Слот палитры расходуется ТОЛЬКО тогда, когда цвет не задан явно —
# ровно то поведение, которое раньше возникало неявно (Makie раздавал цвета по
# порядку создания графиков). Явно заданный цвет (например, синий у заделки или
# оранжевый у ползуна) слота не расходует, поэтому внешний вид сохраняется.
#
# ПРО НЕЗАДАННЫЕ АТРИБУТЫ. `nothing` в `linewidth` / `alpha` / `markersize`
# означает «значение из темы Makie». Такие атрибуты НЕ передаются в Makie вовсе:
# часть графиков (scatter, poly, arrows2d) использует их в арифметике и падает
# на `Makie.automatic`. Поэтому все аргументы собираются через `_add!`, который
# молча пропускает `nothing` и `automatic`.
#
# Атрибут рецепта `Makie.automatic` означает «не переопределять, взять из
# компаньона»: так глобальный `systemfigure(sys, sol; color = :red)` перекрашивает
# всё, а `color = nothing` у отдельной сущности оставляет её на палитре.

# ---------------------------------------------------------------------------
# Раздача цветов из палитры
# ---------------------------------------------------------------------------

"""
    PaletteCycle()

Счётчики палитры на одну ось.

Makie ведёт цикл палитры НЕЗАВИСИМО для каждого семейства графиков: `lines` и
`scatter` в одной оси считаются раздельно (проверено экспериментально), а
`arrows2d` вообще берёт чёрный из темы и палитру не расходует. Поэтому и здесь
счётчик отдельный на семейство.

Внутри композитного рецепта встроенный цикл Makie не работает: все вложенные
графики получают один и тот же цвет. Палитру раздаём сами.
"""
mutable struct PaletteCycle
    colors::Vector
    counters::Dict{Symbol,Int}
end

PaletteCycle() = PaletteCycle(collect(Makie.wong_colors()), Dict{Symbol,Int}())

function _next_palette!(pc::PaletteCycle, family::Symbol)
    i = get(pc.counters, family, 0) + 1
    pc.counters[family] = i
    return pc.colors[mod1(i, length(pc.colors))]
end

"""
    take_color!(pc, family, c)

`c === nothing` — «взять следующий цвет палитры для семейства `family`»;
явно заданный цвет возвращается как есть и слота не расходует (как в Makie).
"""
take_color!(pc::PaletteCycle, family::Symbol, c) =
    c === nothing ? _next_palette!(pc, family) : c

# Makie.automatic («не переопределено») → значение из компаньона
_pick(attr, companion) = attr === Makie.automatic ? companion : attr

# Собрать аргументы, пропустив «не задано» (nothing) и «не переопределено»
# (automatic): такие значения нельзя передавать в Makie явно.
function _add!(kw::Vector{Pair{Symbol,Any}}, key::Symbol, value)
    (value === nothing || value === Makie.automatic) && return kw
    push!(kw, key => value)
    return kw
end

_new_kw() = Pair{Symbol,Any}[]

# Общий набор атрибутов стиля для всех рецептов сущностей
_style_attrs() = (
    time = 1,
    color = Makie.automatic,
    strokecolor = Makie.automatic,
    linewidth = Makie.automatic,
    alpha = Makie.automatic,
    linestyle = Makie.automatic,
)

"""
    _style_kwargs!(kw, p, st) -> kw

Добавить общие атрибуты стиля: атрибут рецепта, а если он не переопределён —
значение из компаньона `st`.
"""
function _style_kwargs!(kw::Vector{Pair{Symbol,Any}}, p, st::VisualStyle)
    _add!(kw, :color, _pick(p.color[], st.color))
    _add!(kw, :linewidth, _pick(p.linewidth[], st.linewidth))
    _add!(kw, :alpha, _pick(p.alpha[], st.alpha))
    _add!(kw, :linestyle, _pick(p.linestyle[], st.linestyle))
    return kw
end

_style_of(v::AbstractVisual) = style_of(v)
_style_of(::NoVisual) = VisualStyle()

# Порядок отрисовки. У пользовательских визуалов поля style может не быть.
_zorder(v::AbstractVisual) = style_of(v).zorder
_zorder(::NoVisual) = 0

# ---------------------------------------------------------------------------
# Тела
# ---------------------------------------------------------------------------

"""
    bodyplot(system, solution, body)

Рецепт тела. Геометрия берётся из компаньона тела ([`BarVisual`](@ref),
[`RectVisual`](@ref), [`CircleVisual`](@ref), [`MeshVisual`](@ref)).

Точка расширения для новой формы тела — метод `_draw_body!(plot, visual)`:
новый тип компаньона + один метод, без правок математики и существующих
рецептов.
"""
Makie.@recipe(BodyPlot, system, solution, body) do scene
    Makie.Attributes(; _style_attrs()...)
end

function Makie.plot!(p::BodyPlot)
    _draw_body!(p, visual(p.body[]))
    return p
end

function _draw_body!(p::BodyPlot, v::AbstractBodyVisual)
    pts = Makie.lift(p.time) do t
        body_contour(p.system[], p.body[], view(p.solution[], :, t))
    end
    kw = _style_kwargs!(_new_kw(), p, v.style)
    Makie.lines!(p, pts; kw...)
    return p
end

# ---------------------------------------------------------------------------
# Соединения
# ---------------------------------------------------------------------------

"""
    hingeplot(system, solution, joint)

Шарнир — маркер в точке соединения.
"""
Makie.@recipe(HingePlot, system, solution, joint) do scene
    Makie.Attributes(; _style_attrs()...,
                     markersize = Makie.automatic,
                     marker = Makie.automatic)
end

function Makie.plot!(p::HingePlot)
    pt = Makie.lift(p.time) do t
        get_hinge_point(p.system[], p.joint[], view(p.solution[], :, t))
    end
    v = visual(p.joint[])
    st = v isa HingeVisual ? v : HingeVisual()
    kw = _new_kw()
    _add!(kw, :color, _pick(p.color[], st.style.color))
    _add!(kw, :markersize, _pick(p.markersize[], st.markersize))
    _add!(kw, :marker, _pick(p.marker[], st.marker))
    _add!(kw, :alpha, _pick(p.alpha[], st.style.alpha))
    Makie.scatter!(p, pt; kw...)
    return p
end

"""
    fixedplot(system, solution, joint)

Заделка — заштрихованная полоса под звеном.
"""
Makie.@recipe(FixedPlot, system, solution, joint) do scene
    Makie.Attributes(; _style_attrs()...)
end

function Makie.plot!(p::FixedPlot)
    pts = Makie.lift(p.time) do t
        fixed_joint_hatch_points(p.system[], p.joint[], view(p.solution[], :, t))
    end
    v = visual(p.joint[])
    st = v isa FixedVisual ? v : FixedVisual()
    kw = _style_kwargs!(_new_kw(), p, st.style)
    Makie.linesegments!(p, pts; kw...)
    return p
end

"""
    sliderplot(system, solution, joint)

Ползун: ось на теле 1 (линия) и прямоугольник на теле 2 (полигон).
Две части стилизуются независимо: `style` — ось, `prism_style` — ползун.
"""
Makie.@recipe(SliderPlot, system, solution, joint) do scene
    Makie.Attributes(; _style_attrs()...,
                     prism_color = Makie.automatic,
                     prism_strokecolor = Makie.automatic)
end

function Makie.plot!(p::SliderPlot)
    axis_pts = Makie.lift(p.time) do t
        slider_axis_points(p.system[], p.joint[], view(p.solution[], :, t))
    end
    prism_pts = Makie.lift(p.time) do t
        slider_prism_points(p.system[], p.joint[], view(p.solution[], :, t))
    end

    v = visual(p.joint[])
    st = v isa SliderVisual ? v : SliderVisual()
    ps = st.prism_style

    # тело 1 — ось соединения (линия направляющей)
    axis_kw = _style_kwargs!(_new_kw(), p, st.style)
    Makie.lines!(p, axis_pts; axis_kw...)

    # тело 2 — ползун
    prism_kw = _new_kw()
    _add!(prism_kw, :color, _pick(p.prism_color[], ps.color))
    _add!(prism_kw, :strokecolor, _pick(p.prism_strokecolor[], ps.strokecolor))
    _add!(prism_kw, :strokewidth, st.prism_strokewidth)
    _add!(prism_kw, :alpha, _pick(p.alpha[], ps.alpha))
    Makie.poly!(p, prism_pts; prism_kw...)
    return p
end

# ---------------------------------------------------------------------------
# Силы
# ---------------------------------------------------------------------------

Makie.@recipe(LinearSpringPlot, system, solution, spring) do scene
    Makie.Attributes(; _style_attrs()...)
end

function Makie.plot!(p::LinearSpringPlot)
    pts = Makie.lift(p.time) do t
        linear_spring_points(p.system[], p.spring[], view(p.solution[], :, t))
    end
    kw = _style_kwargs!(_new_kw(), p, _style_of(visual(p.spring[])))
    Makie.lines!(p, pts; kw...)
    return p
end

Makie.@recipe(TorsionalSpringPlot, system, solution, spring) do scene
    Makie.Attributes(; _style_attrs()...)
end

function Makie.plot!(p::TorsionalSpringPlot)
    pts = Makie.lift(p.time) do t
        torsional_spring_points(p.system[], p.spring[], view(p.solution[], :, t))
    end
    kw = _style_kwargs!(_new_kw(), p, _style_of(visual(p.spring[])))
    Makie.lines!(p, pts; kw...)
    return p
end

"""
    arrowplot(system, solution, force)

Вектор силы — стрелка от точки приложения.
"""
Makie.@recipe(ArrowPlot, system, solution, force) do scene
    Makie.Attributes(; _style_attrs()...,
                     scale = Makie.automatic,
                     head_size = Makie.automatic)
end

function Makie.plot!(p::ArrowPlot)
    pt = Makie.lift(p.time) do t
        first(force_anchor(p.system[], p.force[], view(p.solution[], :, t)))
    end
    dir = Makie.lift(p.time) do t
        last(force_anchor(p.system[], p.force[], view(p.solution[], :, t)))
    end
    v = visual(p.force[])
    st = v isa ArrowVisual ? v : ArrowVisual()
    kw = _new_kw()
    _add!(kw, :color, _pick(p.color[], st.style.color))
    _add!(kw, :alpha, _pick(p.alpha[], st.style.alpha))
    _add!(kw, :lengthscale, _pick(p.scale[], st.scale))
    _add!(kw, :tiplength, _pick(p.head_size[], st.head_size))
    Makie.arrows2d!(p, pt, dir; kw...)
    return p
end

# ---------------------------------------------------------------------------
# Композитный рецепт системы
# ---------------------------------------------------------------------------

"""
    systemplot(system, solution)

Нарисовать всю систему: тела, соединения и силы. Возвращает
`(figure, axis, plot)`; кадр переключается атрибутом `plot.time`.

    using CairoMakie            # или GLMakie, или WGLMakie — рецепт не зависит
    fig, ax, p = systemplot(sys, sol)
    p.time = 500
    save("frame.png", fig)

    record(fig, "out.mp4", 1:10:size(sol, 2)) do t
        p.time = t
    end

Сущности без визуального компаньона (или с `NoVisual`, или со
`style.visible = false`) молча пропускаются — это легальный случай, а не
ошибка: не всё должно рисоваться.

Порядок отрисовки задаётся `style.zorder` компаньона (устойчивая сортировка:
при равных значениях порядок — как при добавлении: тела, соединения, силы).
"""
Makie.@recipe(SystemPlot, system, solution) do scene
    Makie.Attributes(; _style_attrs()...,
                     show_bodies = true,
                     show_connectors = true,
                     show_forces = true)
end

function Makie.plot!(p::SystemPlot)
    sys = p.system[]

    entities = Any[]
    zorders = Int[]
    if p.show_bodies[]
        for bd in bodies(sys)
            is_drawn(bd) || continue
            push!(entities, bd)
            push!(zorders, _zorder(visual(bd)))
        end
    end
    if p.show_connectors[]
        for cn in connectors(sys)
            is_drawn(cn) || continue
            push!(entities, cn)
            push!(zorders, _zorder(visual(cn)))
        end
    end
    if p.show_forces[]
        for fc in sys.forces
            is_drawn(fc) || continue
            push!(entities, fc)
            push!(zorders, _zorder(visual(fc)))
        end
    end

    order = sortperm(zorders; alg = Base.Sort.MergeSort)
    pc = PaletteCycle()
    for i in order
        _plot_entity!(p, pc, entities[i])
    end
    return p
end

# --- выбор рецепта по типу сущности -----------------------------------------

function _plot_entity!(p::SystemPlot, pc::PaletteCycle, bd::Body2D)
    c = take_color!(pc, :lines, _pick(p.color[], _style_of(visual(bd)).color))
    bodyplot!(p, p.system, p.solution, bd;
              time = p.time, color = c, linewidth = p.linewidth,
              alpha = p.alpha, linestyle = p.linestyle)
    return nothing
end

function _plot_entity!(p::SystemPlot, pc::PaletteCycle, joint::HingeJoint)
    c = take_color!(pc, :scatter, _pick(p.color[], _style_of(visual(joint)).color))
    hingeplot!(p, p.system, p.solution, joint;
               time = p.time, color = c, alpha = p.alpha)
    return nothing
end

function _plot_entity!(p::SystemPlot, pc::PaletteCycle, joint::FixedJoint)
    # у заделки цвет по умолчанию задан явно (первый цвет Wong), поэтому слота
    # палитры она не расходует
    c = take_color!(pc, :linesegments, _pick(p.color[], _style_of(visual(joint)).color))
    fixedplot!(p, p.system, p.solution, joint;
               time = p.time, color = c, linewidth = p.linewidth, alpha = p.alpha)
    return nothing
end

function _plot_entity!(p::SystemPlot, pc::PaletteCycle, joint::SliderJoint)
    c = take_color!(pc, :lines, _pick(p.color[], _style_of(visual(joint)).color))
    sliderplot!(p, p.system, p.solution, joint;
                time = p.time, color = c, alpha = p.alpha)
    return nothing
end

function _plot_entity!(p::SystemPlot, pc::PaletteCycle, spring::LinearSpring)
    c = take_color!(pc, :lines, _pick(p.color[], _style_of(visual(spring)).color))
    linearspringplot!(p, p.system, p.solution, spring; time = p.time, color = c)
    return nothing
end

function _plot_entity!(p::SystemPlot, pc::PaletteCycle, spring::TorsionalSpring)
    c = take_color!(pc, :lines, _pick(p.color[], _style_of(visual(spring)).color))
    torsionalspringplot!(p, p.system, p.solution, spring; time = p.time, color = c)
    return nothing
end

function _plot_entity!(p::SystemPlot, pc::PaletteCycle, force::BodyTimeVariableForce)
    # arrows2d! по умолчанию чёрный (значение темы), палитру не расходует:
    # поэтому сюда палитра не подставляется
    c = _pick(p.color[], _style_of(visual(force)).color)
    arrowplot!(p, p.system, p.solution, force; time = p.time, color = c)
    return nothing
end

# Сущность с визуальным компаньоном, для которой рецепта нет, молча
# пропускается: не всё, что существует в модели, должно быть на картинке.
_plot_entity!(p::SystemPlot, pc::PaletteCycle, entity) = nothing

# ---------------------------------------------------------------------------
# Готовая фигура
# ---------------------------------------------------------------------------

"""
    systemfigure(system, solution; limits = :auto, kwargs...) -> (fig, ax, plot)

Создать `Figure` с осью в механическом масштабе (`DataAspect`) и нарисовать
систему. Самый короткий путь:

    using CairoMakie           # или GLMakie, или WGLMakie
    fig, ax, p = systemfigure(sys, sol)
    p.time = 500
    save("frame.png", fig)

Про `limits`:

  * `:auto` (по умолчанию) — границы по всей траектории
    ([`trajectory_limits`](@ref)); кадр стабилен при записи анимации и не нужно
    подбирать числа руками;
  * `nothing` — не трогать границы, оставить Makie;
  * кортеж `(xmin, xmax, ymin, ymax)` — задать явно.

Эквивалент явного варианта, если нужна своя ось:

    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect())
    p = systemplot!(ax, sys, sol)

Масштаб осей нельзя выставить изнутри рецепта (у рецепта нет доступа к объекту
`Axis` — `plot.parent` даёт только `Scene`), поэтому равный масштаб задаётся
здесь.
"""
function systemfigure(sys::MBSystem2D, sol; limits = :auto, kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect())
    p = systemplot!(ax, sys, sol; kwargs...)
    if limits === :auto
        limits!(ax, trajectory_limits(sys, sol)...)
    elseif limits !== nothing
        limits!(ax, limits...)
    end
    return fig, ax, p
end
