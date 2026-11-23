
mutable struct FixedJoint <: AbstractJoint2D
    body::Body2D
    pos::SVector{2,Float64}
    θ::Float64
    index::Int64
    visual::AbstractVisual

    function FixedJoint(body::Body2D; visual::AbstractVisual = FixedVisual())
        return new(body, [0, 0], 0, -1, visual)
    end
end

visual(joint::FixedJoint) = joint.visual

function setposition!(joint::FixedJoint, pos)
    joint.pos = SA[pos[1], pos[2]]
    return nothing
end

function setrotation!(joint::FixedJoint, θ)
    joint.θ = θ
    return nothing
end

number_of_dofs(::FixedJoint) = 3


function get_lms(sys::MBSystem2D, joint::FixedJoint)
    last_lm = sys.lmdofs[joint.index]
    return SA[
        last_lm-2,
        last_lm-1,
        last_lm-0,
    ]
end

function get_fixed_point(sys::MBSystem2D, joint::FixedJoint, state::AbstractVector{Float64})

    bd = joint.body
    pos_dofs = get_body_position_dofs(sys, bd)

    _xi = state[pos_dofs[1]]
    _yi = state[pos_dofs[2]]

    return Point2f(_xi, _yi)
end

function add_to_rhs!(rhs, state, sys::MBSystem2D, joint::FixedJoint)
    body = joint.body
    last_body_dof = sys.bodiesdofs[body.index]

    position_dofs = [
        last_body_dof - 5,
        last_body_dof - 4,
        last_body_dof - 3,
    ]

    velocity_dofs = [
        last_body_dof - 2,
        last_body_dof - 1,
        last_body_dof,
    ]

    joint_dofs = get_lms(sys, joint)

    rhs[velocity_dofs[1]] += state[joint_dofs[1]]
    rhs[velocity_dofs[2]] += state[joint_dofs[2]]
    rhs[velocity_dofs[3]] += state[joint_dofs[3]]

    rhs[joint_dofs[1]] = state[position_dofs[1]] - joint.pos[1]
    rhs[joint_dofs[2]] = state[position_dofs[2]] - joint.pos[2]
    rhs[joint_dofs[3]] = state[position_dofs[3]] - joint.θ
end
function compute_kinematic_residual!(residual::Vector{Float64}, coordinates::Vector{Float64}, sys::MBSystem2D, joint::FixedJoint)
    body = joint.body

    generalized_dofs = get_body_generalized_dofs(sys, body)

    joint_dofs = get_lms(sys, joint) .- last_body_dof(sys)

    residual[joint_dofs[1]] = coordinates[generalized_dofs[1]] - joint.pos[1]
    residual[joint_dofs[2]] = coordinates[generalized_dofs[2]] - joint.pos[2]
    residual[joint_dofs[3]] = coordinates[generalized_dofs[3]] - joint.θ
    return nothing
end

"""
    _eval_drift(f, θ) -> SVector{2}

Смещение точки МЦС шарнира от номинальной точки крепления при относительном угле
`θ`. `nothing` означает «дрейфа нет»: ветвление устраняется на этапе компиляции,
поэтому обычный шарнир не платит за поддержку дрейфа ни аллокациями, ни типовой
нестабильностью.

`θ` не ограничивается типом: при сборке якобиана сюда приходит
`ForwardDiff.Dual`. Результат — двухкомпонентный вектор смещения в ЛСК тела.
"""
@inline _eval_drift(::Nothing, θ) = SA[0.0, 0.0]
@inline _eval_drift(f, θ) = f(θ)

"""
    HingeJoint(bd1, bd2; visual = HingeVisual(), rot_offset = 0.0,
               drift1 = nothing, drift2 = nothing)

Цилиндрический шарнир между телами `bd1` и `bd2`. Номинальные точки крепления
задаются `set_position_on_first_body!` и `set_position_on_second_body!`.

Если заданы `drift1`/`drift2`, точка контакта уходит от номинальной вдоль
аналитических центроид (МЦС): функции получают относительный угол
`θ = θj - θi - rot_offset` и возвращают смещение в ЛСК своего тела — `drift1`
задаёт неподвижную центроиду на `bd1`, `drift2` — подвижную на `bd2`. Начало
отсчёта угла задаёт `rot_offset`, см. `set_rot_offset!`.

Без `drift1`/`drift2` шарнир ведёт себя как обычный: тип
`HingeJoint{Nothing, Nothing}`, смещения нулевые.
"""
mutable struct HingeJoint{F1, F2} <: AbstractJoint2D
    body1::Body2D
    body1_hinge_point::SVector{2,Float64}
    body2::Body2D
    body2_hinge_point::SVector{2,Float64}
    index::Int64
    rot_offset::Float64
    drift1::F1
    drift2::F2
    visual::AbstractVisual

    function HingeJoint(
        bd1::Body2D,
        bd2::Body2D;
        visual::AbstractVisual = HingeVisual(),
        rot_offset::Float64 = 0.0,
        drift1::F1 = nothing,
        drift2::F2 = nothing,
    ) where {F1, F2}
        return new{F1, F2}(
            bd1,
            SA[0.0, 0.0],
            bd2,
            SA[0.0, 0.0],
            -1,
            rot_offset,
            drift1,
            drift2,
            visual,
        )
    end
end

"""
    set_rot_offset!(joint::HingeJoint, offset)

Начало отсчёта относительного угла, который видят `drift1`/`drift2`:
`θ = θj - θi - offset`.
"""
function set_rot_offset!(joint::HingeJoint, offset::Float64)
    joint.rot_offset = offset
    return nothing
end

visual(joint::HingeJoint) = joint.visual

number_of_dofs(::HingeJoint) = 2
function setid!(jnt::T, index) where T<:AbstractJoint2D
    jnt.index = index
end
function set_position_on_first_body!(joint::HingeJoint, pos::SVector{2,Float64})
    joint.body1_hinge_point = pos
    return nothing
end

function set_position_on_second_body!(joint::HingeJoint, pos::SVector{2,Float64})
    joint.body2_hinge_point = pos
    return nothing
end

function get_lms(sys::MBSystem2D, joint::HingeJoint)
    last_lm = sys.lmdofs[joint.index]
    return SA[
        last_lm-1,
        last_lm-0,
    ]
end

function add_to_rhs!(rhs, state, sys::MBSystem2D, joint::HingeJoint)
    bd1 = joint.body1
    bd2 = joint.body2

    bd1_p_dofs = get_body_position_dofs(sys, bd1)
    bd1_v_dofs = get_body_velocity_dofs(sys, bd1)
    bd2_p_dofs = get_body_position_dofs(sys, bd2)
    bd2_v_dofs = get_body_velocity_dofs(sys, bd2)

    _xi = state[bd1_p_dofs[1]]
    _yi = state[bd1_p_dofs[2]]
    _θi = state[bd1_p_dofs[3]]

    _xj = state[bd2_p_dofs[1]]
    _yj = state[bd2_p_dofs[2]]
    _θj = state[bd2_p_dofs[3]]

    lms = get_lms(sys, joint)
    λ1 = state[lms[1]]
    λ2 = state[lms[2]]

    # Точка МЦС в ЛСК каждого тела: номинальная плюс смещение вдоль центроиды.
    # Для обычного шарнира смещение нулевое и ветвление убирается при компиляции.
    θ = _θj - _θi - joint.rot_offset
    d1 = _eval_drift(joint.drift1, θ)
    d2 = _eval_drift(joint.drift2, θ)

    xci = joint.body1_hinge_point[1] + d1[1]
    yci = joint.body1_hinge_point[2] + d1[2]
    xcj = joint.body2_hinge_point[1] + d2[1]
    ycj = joint.body2_hinge_point[2] + d2[2]

    si, ci = sincos(_θi)
    sj, cj = sincos(_θj)

    rxi = xci * ci - yci * si
    ryi = xci * si + yci * ci
    rxj = xcj * cj - ycj * sj
    ryj = xcj * sj + ycj * cj

    rhs[bd1_v_dofs[1]] += λ1
    rhs[bd1_v_dofs[2]] += λ2
    rhs[bd1_v_dofs[3]] += -ryi * λ1 + rxi * λ2

    rhs[bd2_v_dofs[1]] += -λ1
    rhs[bd2_v_dofs[2]] += -λ2
    rhs[bd2_v_dofs[3]] += ryj * λ1 - rxj * λ2

    rhs[lms[1]] = (_xi + rxi) - (_xj + rxj)
    rhs[lms[2]] = (_yi + ryi) - (_yj + ryj)

    return nothing
end


function compute_kinematic_residual!(residual::Vector{Float64}, coordinates::Vector{Float64}, sys::MBSystem2D, joint::HingeJoint)
    bd1 = joint.body1
    bd2 = joint.body2

    bd1_g_dofs = get_body_generalized_dofs(sys, bd1)
    bd2_g_dofs = get_body_generalized_dofs(sys, bd2)

    joint_dofs = get_lms(sys, joint) .- last_body_dof(sys)

    _xi = coordinates[bd1_g_dofs[1]]
    _yi = coordinates[bd1_g_dofs[2]]
    _θi = coordinates[bd1_g_dofs[3]]

    _xj = coordinates[bd2_g_dofs[1]]
    _yj = coordinates[bd2_g_dofs[2]]
    _θj = coordinates[bd2_g_dofs[3]]

    # Точка МЦС берётся так же, как в add_to_rhs!: номинальная плюс дрейф.
    d1 = _eval_drift(joint.drift1, _θj - _θi - joint.rot_offset)
    d2 = _eval_drift(joint.drift2, _θj - _θi - joint.rot_offset)

    xci = joint.body1_hinge_point[1] + d1[1]
    yci = joint.body1_hinge_point[2] + d1[2]
    xcj = joint.body2_hinge_point[1] + d2[1]
    ycj = joint.body2_hinge_point[2] + d2[2]

    residual[joint_dofs[1]] = (_xi + xci * cos(_θi) - yci * sin(_θi)) - 
                         (_xj + xcj * cos(_θj) - ycj * sin(_θj))
    residual[joint_dofs[2]] = (_yi + xci * sin(_θi) + yci * cos(_θi)) -
                         (_yj + xcj * sin(_θj) + ycj * cos(_θj))


    return nothing;
end

function get_hinge_point(system::MBSystem2D, joint::HingeJoint, state::AbstractVector{Float64})
    bd1 = joint.body1
    pos_dofs1 = get_body_position_dofs(system, bd1)
    _xi = state[pos_dofs1[1]]
    _yi = state[pos_dofs1[2]]
    _θi = state[pos_dofs1[3]]

    bd2 = joint.body2
    pos_dofs2 = get_body_position_dofs(system, bd2)
    _θj = state[pos_dofs2[3]]

    # Маркер стоит на фактической точке МЦС, а не на номинальной.
    d1 = _eval_drift(joint.drift1, _θj - _θi - joint.rot_offset)

    xci = joint.body1_hinge_point[1] + d1[1]
    yci = joint.body1_hinge_point[2] + d1[2]

    return Point2f(
        _xi + cos(_θi) * xci - sin(_θi) * yci,
        _yi + sin(_θi) * xci + cos(_θi) * yci
    )
end

mutable struct SliderJoint <: AbstractJoint2D
    body1::Body2D
    body1_position::SVector{2,Float64}
    body1_direction::SVector{2,Float64}
    body2::Body2D
    body2_position::SVector{2,Float64}
    body2_direction::SVector{2,Float64}
    alpha1::Float64
    alpha2::Float64
    index::Int64
    visual::AbstractVisual


    # kwarg'и vis_* сохранены ради совместимости: они конструируют SliderVisual.
    function SliderJoint(bd1::Body2D, bd2::Body2D;
                         vis_axis_length::Float64 = 2.0,
                         vis_slider_length::Float64 = 0.5,
                         vis_slider_width::Float64 = 0.2,
                         visual::AbstractVisual = SliderVisual(axis_length = vis_axis_length,
                                                               slider_length = vis_slider_length,
                                                               slider_width = vis_slider_width))
        α1 = 0;
        α2 = 0;
        return new(bd1, SA[0.0, 0.0], SA[1.0, 0.0], bd2, SA[0.0, 0.0], SA[1.0, 0.0], α1, α2, -1,
                   visual)
    end
end

visual(joint::SliderJoint) = joint.visual

# Параметры отрисовки SliderJoint: совместимость со старым kwarg-API.
# Мутирует вложенный SliderVisual, а не поля сущности.
function set_visualization!(joint::SliderJoint; vis_axis_length::Union{Float64,Nothing} = nothing,
                            vis_slider_length::Union{Float64,Nothing} = nothing,
                            vis_slider_width::Union{Float64,Nothing} = nothing)
    v = visual(joint)
    if (v isa SliderVisual)
        joint.visual = SliderVisual(
            axis_length = vis_axis_length === nothing ? v.axis_length : vis_axis_length,
            slider_length = vis_slider_length === nothing ? v.slider_length : vis_slider_length,
            slider_width = vis_slider_width === nothing ? v.slider_width : vis_slider_width,
            attachment = v.attachment,
            style = v.style,
            prism_style = v.prism_style,
        )
    end
    return nothing
end

number_of_dofs(::SliderJoint) = 2

function get_lms(sys::MBSystem2D, joint::SliderJoint)
    last_lm = sys.lmdofs[joint.index]
    return SA[
        last_lm-1,
        last_lm-0,
    ]
end

function set_position_on_first_body!(joint::SliderJoint, pos::SVector{2,Float64})
    joint.body1_position = pos
    return nothing
end

function set_position_on_second_body!(joint::SliderJoint, pos::SVector{2,Float64})
    joint.body2_position = pos
    return nothing
end

function set_direction_on_first_body!(joint::SliderJoint, dir::SVector{2,Float64})
    joint.alpha1 = atan(dir[2], dir[1])
    joint.body1_direction = dir
    return nothing
end

function set_direction_on_second_body!(joint::SliderJoint, dir::SVector{2,Float64})
    joint.alpha2 = atan(dir[2], dir[1])
    joint.body2_direction = dir
    return nothing
end

function add_to_rhs!(rhs, state, sys::MBSystem2D, joint::SliderJoint)
    bd1 = joint.body1
    bd2 = joint.body2

    bd1_p_dofs = get_body_position_dofs(sys, bd1)
    bd1_v_dofs = get_body_velocity_dofs(sys, bd1)
    bd2_p_dofs = get_body_position_dofs(sys, bd2)
    bd2_v_dofs = get_body_velocity_dofs(sys, bd2)

    _xi = state[bd1_p_dofs[1]]
    _yi = state[bd1_p_dofs[2]]
    _θi = state[bd1_p_dofs[3]]

    # vxi = state[bd1_v_dofs[1]]
    # vyi = state[bd1_v_dofs[2]]
    # ωi = state[bd1_v_dofs[3]]

    _xj = state[bd2_p_dofs[1]]
    _yj = state[bd2_p_dofs[2]]
    _θj = state[bd2_p_dofs[3]]

    # vxj = state[bd2_v_dofs[1]]
    # vyj = state[bd2_v_dofs[2]]
    # ωj = state[bd2_v_dofs[3]]

    lms = get_lms(sys, joint)
    λ1 = state[lms[1]]
    λ2 = state[lms[2]]

    xci = joint.body1_position[1]
    yci = joint.body1_position[2]
    xdi = joint.body1_direction[1]
    ydi = joint.body1_direction[2]

    xcj = joint.body2_position[1]
    ycj = joint.body2_position[2]
    xdj = joint.body2_direction[1]
    ydj = joint.body2_direction[2]

    αi = joint.alpha1
    αj = joint.alpha2
    
    # Constraint equations for slider joint:
    # 1. Perpendicular distance between points is zero
    # 2. Rotation difference (directions aligned)

    # Position of slider points in global coordinates
    xpi = _xi + xci * cos(_θi) - yci * sin(_θi)
    ypi = _yi + xci * sin(_θi) + yci * cos(_θi)

    xpj = _xj + xcj * cos(_θj) - ycj * sin(_θj)
    ypj = _yj + xcj * sin(_θj) + ycj * cos(_θj)


    # normal to direction in terms of first body in global CS
    G_xni = -xdi*sin(_θi) - ydi*cos(_θi)
    G_yni =  xdi*cos(_θi) - ydi*sin(_θi)

    # Constraint 1: perpendicular distance
    rhs[lms[1]] = (xpj - xpi) * G_xni + (ypj - ypi) * G_yni

    # Constraint 2: direction alignment (rotation difference)
    # dx_gj * dy_gi - dy_gj * dx_gi = 0 (cross product)
    rhs[lms[2]] = _θj + αj - (_θi + αi)

    # Velocity constraints
    # Constraint 1 is Φ₁ = (r_pj - r_pi)·n_i, so ∂Φ₁/∂r_i = -n_i and the generalized
    # forces carried by λ₁ are  -λ₁ n_i  on body 1 and  +λ₁ n_i  on body 2.
    # (The signs here used to be the other way round, which made the λ₁ reaction of a
    # slider do spurious virtual work: B = ∂F/∂λ did not match Cᵀ = (∂Φ/∂q)ᵀ, while
    # every other joint in the library did match exactly.)
    rhs[bd1_v_dofs[1]] -= λ1 * G_xni
    rhs[bd1_v_dofs[2]] -= λ1 * G_yni
    rhs[bd1_v_dofs[3]] -= λ1 * ((_yi - ypj) * G_xni + (xpj - _xi) * G_yni)

    rhs[bd2_v_dofs[1]] += λ1 * G_xni
    rhs[bd2_v_dofs[2]] += λ1 * G_yni
    rhs[bd2_v_dofs[3]] += λ1 * ((ypj - _yj ) * G_xni + (_xj - xpj) * G_yni)

    # For constraint 2: d/dt of direction alignment
    rhs[bd1_v_dofs[3]] += -λ2 
    rhs[bd2_v_dofs[3]] += λ2
end

function compute_kinematic_residual!(residual::Vector{Float64}, coordinates::Vector{Float64}, sys::MBSystem2D, joint::SliderJoint)
    bd1 = joint.body1
    bd2 = joint.body2

    bd1_g_dofs = get_body_generalized_dofs(sys, bd1)
    bd2_g_dofs = get_body_generalized_dofs(sys, bd2)

    _xi = coordinates[bd1_g_dofs[1]]
    _yi = coordinates[bd1_g_dofs[2]]
    _θi = coordinates[bd1_g_dofs[3]]

    _xj = coordinates[bd2_g_dofs[1]]
    _yj = coordinates[bd2_g_dofs[2]]
    _θj = coordinates[bd2_g_dofs[3]]


    xci = joint.body1_position[1]
    yci = joint.body1_position[2]
    xdi = joint.body1_direction[1]
    ydi = joint.body1_direction[2]

    xcj = joint.body2_position[1]
    ycj = joint.body2_position[2]

    αi = joint.alpha1
    αj = joint.alpha2
    
    # Position of slider points in global coordinates
    xpi = _xi + xci * cos(_θi) - yci * sin(_θi)
    ypi = _yi + xci * sin(_θi) + yci * cos(_θi)

    xpj = _xj + xcj * cos(_θj) - ycj * sin(_θj)
    ypj = _yj + xcj * sin(_θj) + ycj * cos(_θj)


    # normal to direction in terms of first body in global CS
    G_xni = -xdi*sin(_θi) - ydi*cos(_θi)
    G_yni =  xdi*cos(_θi) - ydi*sin(_θi)


    joint_dofs = get_lms(sys, joint) .- last_body_dof(sys)
    
    # Constraint 1: perpendicular distance
    residual[joint_dofs[1]] = (xpj - xpi) * G_xni + (ypj - ypi) * G_yni

    # Constraint 2: direction alignment (rotation difference)
    # dx_gj * dy_gi - dy_gj * dx_gi = 0 (cross product)
    residual[joint_dofs[2]] = _θj + αj - (_θi + αi)
end