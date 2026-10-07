function Makie.lift(system, solution, body::Body2D, i::Observable)
    return lift(i) do value
        points = Vector{Point2f}(undef, 2)
        points[1:2] .= get_boundary_points(system, body, view(solution, :, value))
        return points
    end
end
function Makie.lift(system, solution, joint::AbstractJoint2D, i::Observable)
end


function Makie.lift(system, solution, joint::HingeJoint, i::Observable)
    p =  lift(i) do value
        point = get_hinge_point(system, joint, view(solution, :, value)) 
        return point;
    end
    return p;
end
# ---------------------------------------------------------------------------
# FixedJoint: заштрихованная область под звеном (обозначение заделки)
# ---------------------------------------------------------------------------

# первый цвет стандартной палитры Makie (Wong) — синий #0072B2
const FIXED_JOINT_COLOR = Makie.wong_colors()[1]

"""
    get_fixed_joint_hatch_points(system, joint::FixedJoint, state) -> Vector{Point2f}

Штриховка заделки: плоский список пар точек (начало, конец) отрезков штриховки.
Полоса расположена под звеном (в локальной системе тела, вдоль локальной оси `-y`),
её длина равна длине звена `body.length`, высота — четверти длины звена.
Штрихи идут под 45° к звену и обрезаны по границам полосы.
"""
function get_fixed_joint_hatch_points(system::MBSystem2D, joint::FixedJoint, state::AbstractVector{Float64})
    bd = joint.body
    pos_dofs = get_body_position_dofs(system, bd)

    _x = state[pos_dofs[1]]
    _y = state[pos_dofs[2]]
    _θ = state[pos_dofs[3]]

    l = bd.length   # длина полосы = длина звена
    h = l / 4       # высота полосы = четверть длины звена

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
    # Шаг выбран так, чтобы штрихи читались как обычная штриховка заделки.
    step = h / 2
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

function Makie.lift(system, solution, joint::FixedJoint, i::Observable)
    return lift(i) do value
        get_fixed_joint_hatch_points(system, joint, view(solution, :, value))
    end
end

# ---------------------------------------------------------------------------
# SliderJoint: ось (линия) на теле 1 + вытянутый вдоль оси прямоугольник на теле 2
# ---------------------------------------------------------------------------

"""
    get_slider_axis(system, joint::SliderJoint, state) -> (Point2f, SVector{2})

Точка соединения на теле 1 и единичное направление оси соединения,
обе величины — в глобальной системе координат.
"""
function get_slider_axis(system::MBSystem2D, joint::SliderJoint, state::AbstractVector{Float64})
    bd1 = joint.body1
    pos_dofs1 = get_body_position_dofs(system, bd1)
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
    get_slider_axis_points(system, joint::SliderJoint, state) -> Vector{Point2f}

Отрезок оси соединения на теле 1 (центр — точка соединения).
"""
function get_slider_axis_points(system::MBSystem2D, joint::SliderJoint, state::AbstractVector{Float64})
    point, d = get_slider_axis(system, joint, state)
    l = joint.vis_axis_length / 2

    return [
        Point2f(point[1] - l * d[1], point[2] - l * d[2]),
        Point2f(point[1] + l * d[1], point[2] + l * d[2]),
    ]
end

"""
    get_slider_prism_points(system, joint::SliderJoint, state) -> Vector{Point2f}

Замкнутый контур прямоугольника на теле 2: вытянут вдоль оси соединения.
"""
function get_slider_prism_points(system::MBSystem2D, joint::SliderJoint, state::AbstractVector{Float64})
    bd2 = joint.body2
    pos_dofs2 = get_body_position_dofs(system, bd2)
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

    l = joint.vis_slider_length / 2
    w = joint.vis_slider_width / 2

    # локальные углы прямоугольника + замыкающая точка
    local_corners = ((-l, -w), (l, -w), (l, w), (-l, w), (-l, -w))

    return [Point2f(xpj + x * c - y * s, ypj + x * s + y * c) for (x, y) in local_corners]
end

function Makie.lift(system, solution, joint::SliderJoint, i::Observable)
    return lift(i) do value
        get_slider_axis_points(system, joint, view(solution, :, value))
    end
end

function lift_slider_prism(system, solution, joint::SliderJoint, i::Observable)
    return lift(i) do value
        get_slider_prism_points(system, joint, view(solution, :, value))
    end
end


function get_torsionalSpring_point(system::MBSystem2D, spring::Union{TorsionalSpring,LinearSpring}, state::AbstractVector{Float64})
    bd1 = spring.hinge.body1
    pos_dofs1 = get_body_position_dofs(system, bd1)
    _xi1 = state[pos_dofs1[1]]
    _yi1 = state[pos_dofs1[2]]
    _θi1 = state[pos_dofs1[3]]

    bd2 = spring.hinge.body2
    pos_dofs2 = get_body_position_dofs(system, bd2)
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

    return Point2f(_xi ,_yi)
end

function get_Spring_point(system::MBSystem2D, spring::LinearSpring, state::AbstractVector{Float64})
    bd1 = spring.joint.body1
    pos_dofs1 = get_body_position_dofs(system, bd1)
    _xi1 = state[pos_dofs1[1]]
    _yi1 = state[pos_dofs1[2]]
    _θi1 = state[pos_dofs1[3]]

    bd2 = spring.joint.body2
    pos_dofs2 = get_body_position_dofs(system, bd2)
    _xi2 = state[pos_dofs2[1]]
    _yi2 = state[pos_dofs2[2]]
    _θi2 = state[pos_dofs2[3]]

    _xi = (_yi1 - _yi2 + tan(_θi2) * _xi2 - tan(_θi1) * _xi1) / (tan(_θi2) - tan(_θi1))
    _yi = tan(_θi2) * (_xi - _xi2) + _yi2

    return Point2f(_xi ,_yi)
end

function Makie.lift(system, solution, spring::LinearSpring, i::Observable)
    p = lift(i) do value
        point = get_Spring_point(system, spring, view(solution, :, value))
        bd1 = spring.joint.body1
        pos_dofs1 = get_body_position_dofs(system, bd1)
        _xi1, _yi1, _θi1 = view(solution, :, value)[pos_dofs1]
        
        _xi1 += spring.joint.body1_position[1]
        _yi1 += spring.joint.body1_position[2]

        bd2 = spring.joint.body2
        pos_dofs2 = get_body_position_dofs(system, bd2)
        _xi2, _yi2, _θi2 = view(solution, :, value)[pos_dofs2]

        _xi2 += spring.joint.body2_position[1]
        _yi2 += spring.joint.body2_position[2]

        t_spring = spring.vis_r / 2
        N = spring.vis_N
        dx = (_xi2 - _xi1) / (N - 2)
        dy = (_yi2 - _yi1) / (N - 2)
        xs = [_xi1; dx/2; ones(N-3)*dx; dx/2]
        ys = [_yi1; dy/2; ones(N-3)*dy; dy/2]

        x_range = cumsum(xs)
        y_range = cumsum(ys)

        points = Vector{Point2f}(undef, N)
        for j in 1:N
            if (j == 1)
                x = _xi1
                y = _yi1
            elseif (j == N)
                x = _xi2
                y = _yi2
            else
                x = x_range[j]
                y = y_range[j] + t_spring * (-1)^j
            end
            points[j] = Point2f(x, y)
        end
        return points;
    end
    return p
end

function Makie.lift(system, solution, spring::TorsionalSpring, i::Observable)
    p = lift(i) do value
        point = get_torsionalSpring_point(system, spring, view(solution, :, value)) 
        bd1 = spring.hinge.body1
        pos_dofs1 = get_body_position_dofs(system, bd1)
        _xi1, _yi1, _θi1 = view(solution, :, value)[pos_dofs1]

        bd2 = spring.hinge.body2
        pos_dofs2 = get_body_position_dofs(system, bd2)
        _xi2, _yi2, _θi2 = view(solution, :, value)[pos_dofs2]

        start_angel = _θi1 + π
        end_angel = _θi2 + 4*π

        r0 = spring.vis_r/2
        r1 = spring.vis_r
        N = 100

        t = LinRange(start_angel, end_angel, N)
        R = LinRange(r0, r1, N)
        x0 = point[1]
        y0 = point[2]

        points = Vector{Point2f}(undef, N)
        for j in 1:N
            x = R[j] * cos(t[j]) + x0  
            y = R[j] * sin(t[j]) + y0
            points[j] = Point2f(x,y)
        end
        return points;
    end
    return p
end

function draw!(ax, joint::AbstractConnector2D, system::MBSystem2D, solution, iter::Observable)
end
function draw!(ax, joint::HingeJoint, system::MBSystem2D, solution, iter::Observable)
    hinge_point = lift(system, solution, joint, iter);
    scatter!(ax, hinge_point);
end

function draw!(ax, joint::FixedJoint, system::MBSystem2D, solution, iter::Observable)
    # заделка — штриховка под звеном; толщина линии как у пружинок (по умолчанию)
    linesegments!(ax, lift(system, solution, joint, iter); color = FIXED_JOINT_COLOR)
end

function draw!(ax, joint::SliderJoint, system::MBSystem2D, solution, iter::Observable)
    # тело 1 — ось соединения (линия направляющей)
    lines!(ax, lift(system, solution, joint, iter))
    # тело 2 — ползун: прямоугольник, вытянутый вдоль оси соединения
    poly!(ax, lift_slider_prism(system, solution, joint, iter);
          color = (:orange, 0.45), strokewidth = 1.0, strokecolor = :orange)
end

function draw!(ax, joint::TorsionalSpring, system::MBSystem2D, solution, iter::Observable)
    hinge_point = lift(system, solution, joint, iter);
    lines!(ax, hinge_point);
end

function draw!(ax, joint::LinearSpring, system::MBSystem2D, solution, iter::Observable)
    hinge_point = lift(system, solution, joint, iter);
    lines!(ax, hinge_point);
end
function draw!(ax, force::AbstractForce2D, system::MBSystem2D, solution, iter::Observable)
end

function draw!(ax, force::BodyTimeVariableForce, system::MBSystem2D, solution, iter::Observable )
    p = lift(iter) do value
        bd = force.body;
        bd_dofs = get_body_position_dofs(system, bd)
        xb, yb, θb = solution[bd_dofs, value]
        rot = [
            cos(θb) -sin(θb);
            sin(θb) cos(θb);
        ]
        point = Point2d([xb, yb] + rot * force.pos) 
        return point
    end

    d = lift(iter) do value
        bd = force.body;
        bd_dofs = get_body_position_dofs(system, bd)

        xb, yb, θb = solution[bd_dofs, value]
        t = solution[end, value]
        rot = [
            cos(θb) -sin(θb);
            sin(θb) cos(θb);
        ]
        
        dir = Point2d(rot * [force.fx(t), force.fy(t)] ) 
        return dir
    end

    arrows2d!(ax, p, d; lengthscale = 0.01)
end

# (отрисовка FixedJoint — штриховка заделки — определена выше, рядом с draw!(HingeJoint))


function animate(sys::MBSystem2D, sol, time_span, filename; framerate=60, limits = (-1, 1, 1, 1))
    fig = Figure()
    iter = Observable(1)
    ax = Axis(fig[1, 1], aspect = DataAspect())

    for body in bodies(sys)
        bar = lift(sys, sol, body, iter)
        lines!(ax, bar)
    end

    for connector in connectors(sys)
        draw!(ax, connector, sys, sol, iter)
    end
    for force in sys.forces
        draw!(ax, force, sys, sol, iter)
    end

    limits!(ax, limits...)

    record(fig, filename, 1:5:length(time_span);
        framerate=framerate / 5) do t
        iter[] = t
    end
end


function draw_static(sys::MBSystem2D, sol;  limits = (-1, 1, -1, 1))
    fig = Figure()
    iter = Observable(1)
    ax = Axis(fig[1, 1], aspect = DataAspect())

    for body in bodies(sys)
        bar = lift(sys, sol, body, iter)
        lines!(ax, bar)
    end

    for connector in connectors(sys)
        draw!(ax, connector, sys, sol, iter)
    end
    for force in sys.forces
        draw!(ax, force, sys, sol, iter)
    end

    limits!(ax, limits...)

    return (fig, ax);
end