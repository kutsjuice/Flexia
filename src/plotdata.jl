# ===========================================================================
# Чистая геометрия отрисовки ("plot data")
# ===========================================================================
#
# Здесь только вычисление точек по состоянию системы — без Makie-графиков,
# без Observable и без бэкенда. Всё это тестируется юнит-тестами.
#
# Геометрические параметры берутся из визуальных компаньонов (`visual(...)`),
# а не из математических типов.

# ---------------------------------------------------------------------------
# Прикрепление фигуры к телу
# ---------------------------------------------------------------------------

"""
    to_global(x, y, θ, att::Attachment, p::Point2f) -> Point2f

Перевести точку `p`, заданную в локальной системе фигуры, в глобальные
координаты тела с учётом смещения и доворота из `att`.
"""
function to_global(x::Real, y::Real, θ::Real, att::Attachment, p::Point2f)
    c, s = cos(θ), sin(θ)
    ox, oy = att.offset[1], att.offset[2]
    gx = x + ox * c - oy * s
    gy = y + ox * s + oy * c
    φ = θ + att.angle
    c2, s2 = cos(φ), sin(φ)
    return Point2f(gx + p[1] * c2 - p[2] * s2,
                   gy + p[1] * s2 + p[2] * c2)
end

# ---------------------------------------------------------------------------
# Тела: контур по визуальному компаньону
# ---------------------------------------------------------------------------

"""
    visual_extent(body) -> Float64

Характерный линейный размер тела (используется, например, для высоты штриховки
заделки, когда она не задана явно).
"""
visual_extent(body::Body2D) = visual_extent(visual(body))
visual_extent(v::BarVisual) = v.length
visual_extent(v::RectVisual) = max(v.width, v.height)
visual_extent(v::CircleVisual) = 2 * v.radius
visual_extent(v::Link2DVisual) = v.length
visual_extent(v::MeshVisual) = _bbox_extent(v.points)
visual_extent(::NoVisual) = 1.0
visual_extent(::AbstractVisual) = 1.0

function _bbox_extent(points::Vector{Point2f})
    isempty(points) && return 1.0
    xs = first.(points)
    ys = last.(points)
    return max(maximum(xs) - minimum(xs), maximum(ys) - minimum(ys))
end

"""
    body_contour(sys, body, state) -> Vector{Point2f}

Контур тела в глобальных координатах — то, что рисуется линией или полигоном.
Диспетчеризация по типу визуального компаньона: новая форма = новый метод.
"""
function body_contour(sys::MBSystem2D, body::Body2D, state::AbstractVector{Float64})
    pos_inds = get_body_position_dofs(sys, body)
    x, y, θ = state[pos_inds]
    return body_contour(x, y, θ, visual(body))
end

# Нет компаньона — нет геометрии
body_contour(x, y, θ, ::NoVisual) = Point2f[]

# Стержень: историческое поведение — отрезок вдоль локальной оси тела
function body_contour(x, y, θ, v::BarVisual)
    l = v.length
    att = v.attachment
    if v.width <= 0
        p1 = to_global(x, y, θ, att, Point2f(-l / 2, 0.0))
        p2 = to_global(x, y, θ, att, Point2f(l / 2, 0.0))
        return Point2f[p1, p2]
    end
    hw = v.width / 2
    return Point2f[
        to_global(x, y, θ, att, Point2f(-l / 2, -hw)),
        to_global(x, y, θ, att, Point2f(l / 2, -hw)),
        to_global(x, y, θ, att, Point2f(l / 2, hw)),
        to_global(x, y, θ, att, Point2f(-l / 2, hw)),
        to_global(x, y, θ, att, Point2f(-l / 2, -hw)),
    ]
end

# Прямоугольник
function body_contour(x, y, θ, v::RectVisual)
    hw, hh = v.width / 2, v.height / 2
    att = v.attachment
    return Point2f[
        to_global(x, y, θ, att, Point2f(-hw, -hh)),
        to_global(x, y, θ, att, Point2f(hw, -hh)),
        to_global(x, y, θ, att, Point2f(hw, hh)),
        to_global(x, y, θ, att, Point2f(-hw, hh)),
        to_global(x, y, θ, att, Point2f(-hw, -hh)),
    ]
end

# Окружность (замкнутый контур)
function body_contour(x, y, θ, v::CircleVisual; n::Int = 32)
    att = v.attachment
    return Point2f[
        to_global(x, y, θ, att, Point2f(v.radius * cos(t), v.radius * sin(t)))
        for t in range(0, 2π; length = n + 1)
    ]
end

# Произвольный контур
function body_contour(x, y, θ, v::MeshVisual)
    isempty(v.points) && return Point2f[]
    att = v.attachment
    pts = Point2f[to_global(x, y, θ, att, p) for p in v.points]
    if v.closed && (pts[1] != pts[end])
        push!(pts, pts[1])
    end
    return pts
end

"""
    body_contour(x, y, θ, v::Link2DVisual; n_arc) -> Vector{Point2f}

Внешняя граница звена со скруглёнными концами — объединения прямоугольника
`length × width` и двух кружков радиуса `width / 2` с центрами в `(±length / 2, 0)`.

Граница обходится по контуру: нижняя грань → правая полуокружность → верхняя
грань → левая полуокружность. Полуокружности касаются граней в углах
прямоугольника, поэтому стыки гладкие (касательная горизонтальна).

Обычно `width < length`; при `width ≥ length` скругления перекрываются и контур
получается вырожденным — геометрия при этом не корректируется молча.
"""
function body_contour(x, y, θ, v::Link2DVisual; n_arc::Int = 16)
    l = v.length
    r = v.width / 2
    att = v.attachment

    local_pts = Point2f[]
    # нижняя грань: слева направо
    push!(local_pts, Point2f(-l / 2, -r))
    push!(local_pts, Point2f(l / 2, -r))
    # правая оконечность: полуокружность вокруг (l/2, 0), от -π/2 до π/2
    # (первая точка совпадает с концом нижней грани — пропускаем)
    for a in range(-π / 2, π / 2; length = n_arc + 1)[2:end]
        push!(local_pts, Point2f(l / 2 + r * cos(a), r * sin(a)))
    end
    # верхняя грань: справа налево
    push!(local_pts, Point2f(-l / 2, r))
    # левая оконечность: полуокружность вокруг (-l/2, 0), от π/2 до 3π/2
    for a in range(π / 2, 3π / 2; length = n_arc + 1)[2:end]
        push!(local_pts, Point2f(-l / 2 + r * cos(a), r * sin(a)))
    end
    # замкнуть контур (последняя точка левой дуги уже совпадает с первой)
    if local_pts[end] != local_pts[1]
        push!(local_pts, local_pts[1])
    end

    return Point2f[to_global(x, y, θ, att, p) for p in local_pts]
end

# ---------------------------------------------------------------------------
# FixedJoint: штриховка заделки (перенесено из visualize.jl без изменений)
# ---------------------------------------------------------------------------

"""
    fixed_joint_hatch_points(sys, joint, state) -> Vector{Point2f}

Штриховка заделки: плоский список пар точек (начало, конец) отрезков штриховки.
Полоса расположена под звеном (в локальной системе тела, вдоль локальной оси `-y`),
её длина равна длине звена, высота — четверти длины звена.
Штрихи идут под 45° к звену и обрезаны по границам полосы.
"""
function fixed_joint_hatch_points(sys::MBSystem2D, joint::FixedJoint, state::AbstractVector{Float64})
    bd = joint.body
    v = visual(joint)
    pos_dofs = get_body_position_dofs(sys, bd)

    _x = state[pos_dofs[1]]
    _y = state[pos_dofs[2]]
    _θ = state[pos_dofs[3]]

    l = v.hatch_length > 0 ? v.hatch_length : visual_extent(bd)  # длина полосы = длина звена
    h = l * v.hatch_height_ratio                                # высота полосы

    points = Point2f[]
    if (l <= 0 || h <= 0)
        return points
    end

    cθ = cos(_θ)
    sθ = sin(_θ)

    # локальные координаты (s, t): s — вдоль звена, t — вниз, под звеном
    to_global = (s, t) -> Point2f(_x + s * cθ + t * sθ,
                                  _y + s * sθ - t * cθ)

    # Штрих под 45° к звену: в локальных координатах s + t = c.
    # Полоса: s ∈ [-l/2, l/2], t ∈ [0, h] — отсюда s ∈ [c-h, c].
    step = h * v.hatch_step_ratio
    c = -l / 2
    while (c <= l / 2 + h)
        s_lo = max(-l / 2, c - h)
        s_hi = min(l / 2, c)
        if (s_hi > s_lo)
            push!(points, to_global(s_lo, c - s_lo))
            push!(points, to_global(s_hi, c - s_hi))
        end
        c += step
    end

    return points
end

# ---------------------------------------------------------------------------
# SliderJoint: ось на теле 1 + вытянутый прямоугольник на теле 2
# ---------------------------------------------------------------------------

"""
    slider_axis(sys, joint, state) -> (Point2f, SVector{2})

Точка соединения на теле 1 и единичное направление оси соединения,
обе величины — в глобальной системе координат.
"""
function slider_axis(sys::MBSystem2D, joint::SliderJoint, state::AbstractVector{Float64})
    bd1 = joint.body1
    pos_dofs1 = get_body_position_dofs(sys, bd1)
    _xi1 = state[pos_dofs1[1]]
    _yi1 = state[pos_dofs1[2]]
    _θi1 = state[pos_dofs1[3]]

    xci = joint.body1_position[1]
    yci = joint.body1_position[2]

    # точка соединения на теле 1 в глобальных координатах
    xpi = _xi1 + xci * cos(_θi1) - yci * sin(_θi1)
    ypi = _yi1 + xci * sin(_θi1) + yci * cos(_θi1)

    # ось соединения жёстко связана с телом 1 (угол alpha1 задан локально)
    _φi = _θi1 + joint.alpha1
    d = SA_F64[cos(_φi), sin(_φi)]

    return Point2f(xpi, ypi), d
end

"""
    slider_axis_points(sys, joint, state) -> Vector{Point2f}

Отрезок оси соединения на теле 1 (центр — точка соединения).
"""
function slider_axis_points(sys::MBSystem2D, joint::SliderJoint, state::AbstractVector{Float64})
    point, d = slider_axis(sys, joint, state)
    l = visual(joint).axis_length / 2

    return [
        Point2f(point[1] - l * d[1], point[2] - l * d[2]),
        Point2f(point[1] + l * d[1], point[2] + l * d[2]),
    ]
end

"""
    slider_prism_points(sys, joint, state) -> Vector{Point2f}

Замкнутый контур прямоугольника на теле 2: вытянут вдоль оси соединения.
"""
function slider_prism_points(sys::MBSystem2D, joint::SliderJoint, state::AbstractVector{Float64})
    bd2 = joint.body2
    pos_dofs2 = get_body_position_dofs(sys, bd2)
    _xj = state[pos_dofs2[1]]
    _yj = state[pos_dofs2[2]]
    _θj = state[pos_dofs2[3]]

    xcj = joint.body2_position[1]
    ycj = joint.body2_position[2]

    # точка соединения на теле 2 в глобальных координатах
    xpj = _xj + xcj * cos(_θj) - ycj * sin(_θj)
    ypj = _yj + xcj * sin(_θj) + ycj * cos(_θj)

    # ось соединения в системе тела 2
    _φj = _θj + joint.alpha2
    c = cos(_φj)
    s = sin(_φj)

    v = visual(joint)
    l = v.slider_length / 2
    w = v.slider_width / 2

    # локальные углы прямоугольника + замыкающая точка
    local_corners = ((-l, -w), (l, -w), (l, w), (-l, w), (-l, -w))

    return [Point2f(xpj + x * c - y * s, ypj + x * s + y * c) for (x, y) in local_corners]
end

# ---------------------------------------------------------------------------
# Силы: точки приложения и геометрия пружин
# ---------------------------------------------------------------------------

function torsional_spring_anchor(sys::MBSystem2D, spring::Union{TorsionalSpring,LinearSpring}, state::AbstractVector{Float64})
    bd1 = spring.hinge.body1
    pos_dofs1 = get_body_position_dofs(sys, bd1)
    _xi1 = state[pos_dofs1[1]]
    _yi1 = state[pos_dofs1[2]]
    _θi1 = state[pos_dofs1[3]]

    bd2 = spring.hinge.body2
    pos_dofs2 = get_body_position_dofs(sys, bd2)
    _xi2 = state[pos_dofs2[1]]
    _yi2 = state[pos_dofs2[2]]
    _θi2 = state[pos_dofs2[3]]
    R = SA_F64[
        cos(_θi1) -sin(_θi1);
        sin(_θi1) cos(_θi1)
    ]
    v = R * spring.hinge.body1_hinge_point
    _xi = _xi1 + v[1]
    _yi = _yi1 + v[2]

    return Point2f(_xi, _yi)
end

function linear_spring_contact_point(sys::MBSystem2D, spring::LinearSpring, state::AbstractVector{Float64})
    bd1 = spring.joint.body1
    pos_dofs1 = get_body_position_dofs(sys, bd1)
    _xi1 = state[pos_dofs1[1]]
    _yi1 = state[pos_dofs1[2]]
    _θi1 = state[pos_dofs1[3]]

    bd2 = spring.joint.body2
    pos_dofs2 = get_body_position_dofs(sys, bd2)
    _xi2 = state[pos_dofs2[1]]
    _yi2 = state[pos_dofs2[2]]
    _θi2 = state[pos_dofs2[3]]

    _xi = (_yi1 - _yi2 + tan(_θi2) * _xi2 - tan(_θi1) * _xi1) / (tan(_θi2) - tan(_θi1))
    _yi = tan(_θi2) * (_xi - _xi2) + _yi2

    return Point2f(_xi, _yi)
end

"""
    linear_spring_points(sys, spring, state) -> Vector{Point2f}

Зигзаг линейной пружины между точками стыка ползуна.

Концы пружины — точки `body1_position` / `body2_position`, пересчитанные в
глобальную систему **с учётом поворота тел** (раньше локальные координаты
прибавлялись без поворота, и при повёрнутом ползуне пружина рисовалась не на
месте). Зигзаг откладывается по нормали к оси пружины, амплитуда
`coil_radius / 2`, поэтому полная высота зигзага равна `coil_radius` и её можно
подобрать под толщину звена.
"""
function linear_spring_points(sys::MBSystem2D, spring::LinearSpring, state::AbstractVector{Float64})
    j = spring.joint

    contact_point(body, local_pos) = begin
        pos_dofs = get_body_position_dofs(sys, body)
        _x, _y, _θ = state[pos_dofs]
        (x = _x + local_pos[1] * cos(_θ) - local_pos[2] * sin(_θ),
         y = _y + local_pos[1] * sin(_θ) + local_pos[2] * cos(_θ))
    end

    p1 = contact_point(j.body1, j.body1_position)
    p2 = contact_point(j.body2, j.body2_position)

    v = visual(spring)
    amplitude = v.coil_radius / 2      # полная высота зигзага = coil_radius
    N = max(v.n_coils, 4)

    dx = p2.x - p1.x
    dy = p2.y - p1.y
    len = hypot(dx, dy)
    if len < eps(Float64) || amplitude <= 0
        return [Point2f(p1.x, p1.y), Point2f(p2.x, p2.y)]
    end

    ux, uy = dx / len, dy / len
    nx, ny = -uy, ux                    # нормаль к оси пружины

    # Осевой шаг между соседними точками зигзага. Концевые участки идут от оси
    # к первому «гребню», то есть по ширине проходят только половину амплитуды,
    # поэтому по длине им нужно полшага — иначе их наклон вдвое меньше, чем у
    # внутренних участков, и пружина выглядит «с отогнутыми концами».
    Δs = 1.0 / (N - 2)

    points = Vector{Point2f}(undef, N)
    points[1] = Point2f(p1.x, p1.y)
    points[N] = Point2f(p2.x, p2.y)
    for k in 2:(N - 1)
        s = Δs / 2 + (k - 2) * Δs       # полшага до первого гребня, далее полный шаг
        sign = isodd(k) ? 1.0 : -1.0
        points[k] = Point2f(p1.x + s * dx + sign * amplitude * nx,
                            p1.y + s * dy + sign * amplitude * ny)
    end
    return points
end

"""
    torsional_spring_points(sys, spring, state) -> Vector{Point2f}

Спираль крутильной пружины вокруг шарнира.
"""
function torsional_spring_points(sys::MBSystem2D, spring::TorsionalSpring, state::AbstractVector{Float64})
    point = torsional_spring_anchor(sys, spring, state)
    bd1 = spring.hinge.body1
    pos_dofs1 = get_body_position_dofs(sys, bd1)
    _xi1, _yi1, _θi1 = state[pos_dofs1]

    bd2 = spring.hinge.body2
    pos_dofs2 = get_body_position_dofs(sys, bd2)
    _xi2, _yi2, _θi2 = state[pos_dofs2]

    start_angel = _θi1 + π
    end_angel = _θi2 + 4 * π

    v = visual(spring)
    r0 = v.coil_radius / 2
    r1 = v.coil_radius
    N = v.n_points

    t = LinRange(start_angel, end_angel, N)
    R = LinRange(r0, r1, N)
    x0 = point[1]
    y0 = point[2]

    points = Vector{Point2f}(undef, N)
    for j in 1:N
        x = R[j] * cos(t[j]) + x0
        y = R[j] * sin(t[j]) + y0
        points[j] = Point2f(x, y)
    end
    return points
end

"""
    force_anchor(sys, force::BodyTimeVariableForce, state) -> (Point2f, Point2f)

Точка приложения и вектор силы в глобальных координатах.
"""
function force_anchor(sys::MBSystem2D, force::BodyTimeVariableForce, state::AbstractVector{Float64})
    bd = force.body
    bd_dofs = get_body_position_dofs(sys, bd)
    xb, yb, θb = state[bd_dofs]
    rot = [
        cos(θb) -sin(θb);
        sin(θb) cos(θb);
    ]
    point = Point2d([xb, yb] + rot * force.pos)
    t = state[end]
    dir = Point2d(rot * [force.fx(t), force.fy(t)])
    return point, dir
end

# ---------------------------------------------------------------------------
# Автолимиты
# ---------------------------------------------------------------------------

"""
    trajectory_limits(sys, sol; margin = 0.05, stride = 1)
        -> (xmin, xmax, ymin, ymax)

Габарит всей траектории: bounding box по геометрии всех рисуемых тел во всех
кадрах решения плюс отступ `margin` (доля от размера).

Зачем не `autolimits!(ax)`: тот пересчитывает границы по текущему содержимому
на каждом кадре, поэтому при записи анимации картинка дрожит и зумится.
Здесь границы считаются один раз по всему решению и остаются стабильными.

Считаются только тела — стрелки сил могут выходить за отступ.
"""
function trajectory_limits(sys::MBSystem2D, sol::AbstractMatrix; margin = 0.05, stride = 1)
    xs = Float64[]
    ys = Float64[]
    for t in 1:stride:size(sol, 2)
        for bd in bodies(sys)
            is_drawn(bd) || continue
            for pt in body_contour(sys, bd, view(sol, :, t))
                push!(xs, pt[1])
                push!(ys, pt[2])
            end
        end
    end
    isempty(xs) && return (-1.0, 1.0, -1.0, 1.0)

    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)
    # вырожденный случай: всё стоит в одной точке
    if xmax - xmin <= eps(Float64)
        xmin -= 0.5
        xmax += 0.5
    end
    if ymax - ymin <= eps(Float64)
        ymin -= 0.5
        ymax += 0.5
    end
    dx = (xmax - xmin) * margin
    dy = (ymax - ymin) * margin
    return (xmin - dx, xmax + dx, ymin - dy, ymax + dy)
end
