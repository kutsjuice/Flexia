# ===========================================================================
# Визуальные компаньоны (visual companions)
# ===========================================================================
#
# Компаньон — это отдельный объект, описывающий, КАК выглядит физическая
# сущность (тело, соединение, сила). Математические типы (Body2D, HingeJoint,
# LinearSpring, ...) при этом не несут графических полей: единственное, что их
# связывает — ссылка `visual` на компаньон.
#
# Зачем:
#   * тело — не всегда «палка»: BarVisual, RectVisual, CircleVisual, MeshVisual;
#   * новые формы добавляются без правок математики и существующих рецептов;
#   * графические параметры больше не размазаны по математическим типам
#     (раньше: Body2D.length, SliderJoint.vis_axis_length, LinearSpring.vis_r, ...).
#
# Компаньон описывает ТОЛЬКО то, что не меняется во времени. Всё, что зависит от
# кадра (x, y, θ тела), берётся из решения `solution` по атрибуту `time` рецепта.

# ---------------------------------------------------------------------------
# Базовые типы
# ---------------------------------------------------------------------------

"""
    AbstractVisual

Супертип всех визуальных компаньонов.
"""
abstract type AbstractVisual end

"""
    NoVisual()

Сущность рисовать не нужно — это легальное намерение, а не отсутствие метода.

Именно этот тип отличает «датчик/актуатор, который по замыслу не рисуется» от
«забыли реализовать отрисовку»: первое молча пропускается, второе — ошибка.
"""
struct NoVisual <: AbstractVisual end

abstract type AbstractBodyVisual <: AbstractVisual end
abstract type AbstractConnectorVisual <: AbstractVisual end
abstract type AbstractForceVisual <: AbstractVisual end

# ---------------------------------------------------------------------------
# Общие кирпичи: прикрепление и стиль
# ---------------------------------------------------------------------------

"""
    Attachment(; offset, angle)

Где и как фигура прикреплена к телу. Координаты — в локальной системе тела.

  * `offset` — смещение фигуры относительно начала локальной системы тела
    (например, «палка прикреплена не центром, а левым концом»);
  * `angle`  — локальный доворот фигуры относительно тела.
"""
Base.@kwdef struct Attachment
    offset::SVector{2,Float64} = SA[0.0, 0.0]
    angle::Float64 = 0.0
end

"""
    VisualStyle(; color, strokecolor, linewidth, alpha, linestyle,
                  visible, zorder, label)

Стиль отрисовки компаньона.

`color`, `strokecolor`, `linestyle` намеренно нетипизированы (`Any`): Makie
принимает много разных представлений цвета.

**`color = nothing` означает «взять цвет из палитры Makie»** — то самое
поведение, которое было раньше неявным: Makie раздавал цвета по порядку
создания графиков. Явно заданный цвет слота палитры не расходует.
`linewidth`/`alpha`/`markersize` со значением `nothing` означают «значение
по умолчанию из темы Makie».

Любое из этих полей может быть не константой, а функцией от состояния —
см. [`resolve_style`](@ref). Это заложено для раскраски по напряжениям,
скорости и т.п.; поведение по умолчанию — обычный цвет.

  * `visible` — спрятать конкретный экземпляр, не меняя его тип;
  * `zorder`  — порядок наложения. В 2D Makie порядок задаётся порядком
    добавления графиков в ось, поэтому `zorder` реализуется сортировкой
    сущностей внутри рецепта системы;
  * `label`   — подпись для легенды/аннотации.
"""
Base.@kwdef struct VisualStyle
    color::Any = nothing
    strokecolor::Any = nothing
    linewidth::Union{Nothing,Float64} = nothing
    alpha::Union{Nothing,Float64} = nothing
    linestyle::Any = nothing
    visible::Bool = true
    zorder::Int = 0
    label::Union{Nothing,String} = nothing
end

"""
    resolve_style(value, sys, sol, t, entity)

Разрешить значение атрибута стиля: константа возвращается как есть, функция —
вызывается как `value(sys, sol, t, entity)`.

Это единственная точка, через которую рецепты получают стиль, поэтому
поддержка «цвета по состоянию» включается здесь, не трогая рецепты.
"""
resolve_style(value, sys, sol, t, entity) = value
resolve_style(f::Function, sys, sol, t, entity) = f(sys, sol, t, entity)

# ---------------------------------------------------------------------------
# Компаньоны тел
# ---------------------------------------------------------------------------

"""
    BarVisual(; length, width, attachment, style)

Стержень — текущее (историческое) представление тела: отрезок длины `length`,
построенный вдоль локальной оси тела и отцентрованный в начале координат.

`width == 0` — рисуем линией; `width > 0` — прямоугольником такой ширины.
"""
Base.@kwdef struct BarVisual <: AbstractBodyVisual
    length::Float64 = 1.0
    width::Float64 = 0.0
    attachment::Attachment = Attachment()
    style::VisualStyle = VisualStyle()
end

"""
    RectVisual(; width, height, attachment, style)

Прямоугольник `width × height` в локальной системе тела.
"""
Base.@kwdef struct RectVisual <: AbstractBodyVisual
    width::Float64 = 0.2
    height::Float64 = 0.2
    attachment::Attachment = Attachment()
    style::VisualStyle = VisualStyle()
end

"""
    CircleVisual(; radius, attachment, style)

Окружность радиуса `radius` с центром в точке прикрепления.
"""
Base.@kwdef struct CircleVisual <: AbstractBodyVisual
    radius::Float64 = 0.1
    attachment::Attachment = Attachment()
    style::VisualStyle = VisualStyle()
end

"""
    MeshVisual(; points, closed, filled, attachment, style)

Произвольный 2D-контур, заданный в локальной системе тела.

`closed` — замыкать ли контур (добавлять первую точку в конец);
`filled` — заливать (`poly!`) или рисовать только контур (`lines!`).

Этот же тип — точка расширения для раскраски по напряжениям: цвет контура или
заливки можно задать функцией через `style.color`.
"""
Base.@kwdef struct MeshVisual <: AbstractBodyVisual
    points::Vector{Point2f} = Point2f[]
    closed::Bool = true
    filled::Bool = true
    attachment::Attachment = Attachment()
    style::VisualStyle = VisualStyle()
end

"""
    Link2DVisual(; length, width, attachment, style)

Звено с скруглёнными концами: булево объединение прямоугольника
`length × width` и двух кружков радиуса `width / 2`, центры которых лежат
строго на серединах торцевых граней прямоугольника, то есть в
`(±length / 2, 0)`. Диаметр кружка равен длине торцевой грани (`width`).

Центры скруглений совпадают с концами стержня [`BarVisual`](@ref) той же длины —
с типовыми точками крепления шарниров. Поэтому `length` задаёт расстояние между
центрами скруглений (между пальцами), а габарит по длине равен `length + width`.

Имя `Link2DVisual`, а не `2DLink`: идентификатор в Julia не может начинаться с цифры.
"""
Base.@kwdef struct Link2DVisual <: AbstractBodyVisual
    length::Float64 = 1.0
    width::Float64 = 0.2
    attachment::Attachment = Attachment()
    style::VisualStyle = VisualStyle()
end

# ---------------------------------------------------------------------------
# Компаньоны соединений
# ---------------------------------------------------------------------------

"""
    HingeVisual(; markersize, marker, attachment, style)

Шарнир — маркер в точке соединения.
"""
Base.@kwdef struct HingeVisual <: AbstractConnectorVisual
    markersize::Union{Nothing,Float64} = nothing   # nothing → как в теме Makie
    marker::Any = :circle
    attachment::Attachment = Attachment()
    style::VisualStyle = VisualStyle()             # цвет из палитры, как раньше
end

"""
    FixedVisual(; hatch_length, hatch_height_ratio, hatch_angle,
                  hatch_step_ratio, attachment, style)

Заделка — заштрихованная полоса под звеном.

`hatch_length == 0` — брать длину звена (историческое поведение).
Штрихи идут под 45° — угол пока не параметризован.
"""
Base.@kwdef struct FixedVisual <: AbstractConnectorVisual
    hatch_length::Float64 = 0.0
    hatch_height_ratio::Float64 = 0.25
    hatch_step_ratio::Float64 = 0.5
    attachment::Attachment = Attachment()
    # исторический FIXED_JOINT_COLOR = первый цвет палитры Wong
    style::VisualStyle = VisualStyle(color = Makie.wong_colors()[1])
end

"""
    SliderVisual(; axis_length, slider_length, slider_width,
                   style, prism_style)

Ползун: ось (отрезок) на теле 1 плюс вытянутый прямоугольник на теле 2.

Поля — прямые наследники `SliderJoint.vis_axis_length`,
`vis_slider_length`, `vis_slider_width`.
"""
Base.@kwdef struct SliderVisual <: AbstractConnectorVisual
    axis_length::Float64 = 2.0
    slider_length::Float64 = 0.5
    slider_width::Float64 = 0.2
    attachment::Attachment = Attachment()
    style::VisualStyle = VisualStyle()             # ось: цвет из палитры, как раньше
    prism_style::VisualStyle = VisualStyle(color = (:orange, 0.45), strokecolor = :orange)
    prism_strokewidth::Union{Nothing,Float64} = 1.0
end

# ---------------------------------------------------------------------------
# Компаньоны сил
# ---------------------------------------------------------------------------

"""
    LinearSpringVisual(; coil_radius, n_coils, attachment, style)

Линейная пружина: зигзаг вдоль оси ползуна (`vis_r`, `vis_N` исторически).
"""
Base.@kwdef struct LinearSpringVisual <: AbstractForceVisual
    coil_radius::Float64 = 0.3
    n_coils::Int = 4
    attachment::Attachment = Attachment()
    style::VisualStyle = VisualStyle()
end

"""
    TorsionalSpringVisual(; coil_radius, n_points, attachment, style)

Крутильная пружина: спираль вокруг шарнира (исторический `vis_r`).
"""
Base.@kwdef struct TorsionalSpringVisual <: AbstractForceVisual
    coil_radius::Float64 = 0.3
    n_points::Int = 100
    attachment::Attachment = Attachment()
    style::VisualStyle = VisualStyle()
end

"""
    ArrowVisual(; scale, head_size, attachment, style)

Вектор силы — стрелка от точки приложения.
"""
Base.@kwdef struct ArrowVisual <: AbstractForceVisual
    scale::Float64 = 0.01
    head_size::Union{Nothing,Float64} = nothing    # nothing → как в теме Makie
    attachment::Attachment = Attachment()
    style::VisualStyle = VisualStyle()             # цвет из палитры, как раньше
end

# ---------------------------------------------------------------------------
# Доступ к компаньону сущности
# ---------------------------------------------------------------------------

"""
    visual(entity) -> AbstractVisual

Компаньон сущности. По умолчанию — [`NoVisual`](@ref), то есть «не рисовать»:
датчики, актуаторы и прочие нерисуемые сущности не требуют никаких методов
отрисовки и молча пропускаются.
"""
visual(::Any) = NoVisual()

"""
    set_visualization!(entity, v::AbstractVisual) -> nothing

Назначить сущности визуальный компаньон.
"""
function set_visualization!(entity, v::AbstractVisual)
    entity.visual = v
    return nothing
end

"""
    is_drawn(entity) -> Bool

Нужно ли рисовать сущность (есть компаньон, отличный от `NoVisual`, и он видим).

Компаньон может не иметь поля `style` (пользовательский тип с собственной
отрисовкой) — тогда он считается видимым.
"""
is_drawn(entity) = is_drawn_visual(visual(entity))
is_drawn_visual(::NoVisual) = false
is_drawn_visual(v::AbstractVisual) = style_of(v).visible

"""
    style_of(v::AbstractVisual) -> VisualStyle

Стиль компаньона. Если у пользовательского визуала нет поля `style`, берётся
стиль по умолчанию — так свои формы тел не обязаны повторять структуру
встроенных.
"""
style_of(v::AbstractVisual) = hasproperty(v, :style) ? v.style : VisualStyle()
