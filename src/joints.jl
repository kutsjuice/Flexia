
mutable struct FixedJoint <: AbstractJoint2D
    body::Body2D
    pos::SVector{2,Float64}
    θ::Float64
    index::Int64

    function FixedJoint(body::Body2D)
        return new(body, [0, 0], 0, -1)
    end
end

function set_position!(joint::FixedJoint, pos)
    joint.pos = SA[pos[1], pos[2]]
    return nothing
end

function set_rotation!(joint::FixedJoint, θ)
    joint.θ = θ
    return nothing
end

number_of_dofs(::FixedJoint) = 3

function add!(sys::MBSystem2D, joint::AbstractJoint2D)
    push!(sys.joints, joint)
    # sys.jointsnum += 1;
    last_joint_dof = last_lm_dof(sys) + number_of_dofs(joint)
    push!(sys.lmdofs, last_joint_dof)
    setid!(joint, length(sys.joints))
end

function get_lms(sys::MBSystem2D, joint::FixedJoint)
    last_lm = sys.lmdofs[joint.index]
    return SA[
        last_lm-2,
        last_lm-1,
        last_lm-0,
    ]
end

function add_joint_to_rhs!(rhs, state, sys::MBSystem2D, joint::FixedJoint)
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


mutable struct HingeJoint <: AbstractJoint2D
    body1::Body2D
    body1_hinge_point::SVector{2,Float64}
    body2::Body2D
    body2_hinge_point::SVector{2,Float64}
    index::Int64

    function HingeJoint(bd1::Body2D, bd2::Body2D)
        return new(bd1, SA[0.0, 0.0], bd2, SA[0.0, 0.0])
    end
end

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

function add_joint_to_rhs!(rhs, state, sys::MBSystem2D, joint::HingeJoint)
    bd1 = joint.body1
    bd2 = joint.body2

    bd1_p_dofs = get_body_position_dofs(sys, bd1)
    bd1_v_dofs = get_body_velocity_dofs(sys, bd1)
    bd2_p_dofs = get_body_position_dofs(sys, bd2)
    bd2_v_dofs = get_body_velocity_dofs(sys, bd2)

    _xi = state[bd1_p_dofs[1]]
    _yi = state[bd1_p_dofs[2]]
    _θi = state[bd1_p_dofs[3]]

    vxi = state[bd1_v_dofs[1]]
    vyi = state[bd1_v_dofs[2]]
    _ωi = state[bd1_v_dofs[3]]

    _xj = state[bd2_p_dofs[1]]
    _yj = state[bd2_p_dofs[2]]
    _θj = state[bd2_p_dofs[3]]

    vxj = state[bd2_v_dofs[1]]
    vyj = state[bd2_v_dofs[2]]
    _ωj = state[bd2_v_dofs[3]]

    lms = get_lms(sys, joint)
    λ1 = state[lms[1]]
    λ2 = state[lms[2]]

    xci = joint.body1_hinge_point[1]
    yci = joint.body1_hinge_point[2]
    xcj = joint.body2_hinge_point[1]
    ycj = joint.body2_hinge_point[2]

    rhs[bd1_v_dofs[1]] += λ1
    rhs[bd1_v_dofs[2]] += λ2
    rhs[bd1_v_dofs[3]] += (-xci * sin(_θi) - yci * cos(_θi)) * λ1 +
                          (xci * cos(_θi) - yci * sin(_θi)) * λ2

    rhs[bd2_v_dofs[1]] += -λ1
    rhs[bd2_v_dofs[2]] += -λ2
    rhs[bd2_v_dofs[3]] += (xcj * sin(_θj) + ycj * cos(_θj)) * λ1 +
                          (-xcj * cos(_θj) + ycj * sin(_θj)) * λ2

    rhs[lms[1]] = (_xi + xci * cos(_θi) - yci * sin(_θi)) -
                  (_xj + xcj * cos(_θj) - ycj * sin(_θj))
    rhs[lms[2]] = (_yi + xci * sin(_θi) + yci * cos(_θi)) -
                  (_yj + xcj * sin(_θj) + ycj * cos(_θj))


end

function get_hinge_point(system::MBSystem2D, joint::HingeJoint, state::AbstractVector{Float64})
    bd = joint.body1
    pos_dofs = get_body_position_dofs(system, bd)
    _xi = state[pos_dofs[1]]
    _yi = state[pos_dofs[2]]
    _θi = state[pos_dofs[3]]

    xci = joint.body1_hinge_point[1]
    yci = joint.body1_hinge_point[2]

    return Point2f(
        _xi + cos(_θi) * xci - sin(_θi) * yci,
        _yi + sin(_θi) * xci + cos(_θi) * yci
    )
end

mutable struct TorsionalSpring <: AbstractJoint2D
    body1::Body2D
    body2::Body2D
    stiffness::Float64
    rest_angle::Float64
    damping::Float64
    index::Int64
    
    function TorsionalSpring(body1::Body2D, body2::Body2D, 
                            stiffness=1.0, rest_angle=0.0, damping=0.0)
        return new(body1, body2, stiffness, rest_angle, damping, -1)
    end
end

number_of_dofs(::TorsionalSpring) = 0  # пружина не добавляет лагранжевых множителей

function add!(sys::MBSystem2D, spring::AbstractJoint2D)
    push!(sys.joints, spring)
    setid!(spring, length(sys.joints))
end

function add_spring_to_rhs!(rhs, state, sys::MBSystem2D, spring::TorsionalSpring)
    bd1 = spring.body1
    bd2 = spring.body2
    
    # Получаем индексы DOF
    bd1_pos_dofs = get_body_position_dofs(sys, bd1)
    bd1_vel_dofs = get_body_velocity_dofs(sys, bd1)
    bd2_pos_dofs = get_body_position_dofs(sys, bd2)
    bd2_vel_dofs = get_body_velocity_dofs(sys, bd2)
    
    # Текущие углы и угловые скорости
    θ1 = state[bd1_pos_dofs[3]]
    θ2 = state[bd2_pos_dofs[3]]
    ω1 = state[bd1_vel_dofs[3]]
    ω2 = state[bd2_vel_dofs[3]]
    
    # Относительное смещение и скорость
    Δθ = θ1 - θ2 - spring.rest_angle
    Δω = ω1 - ω2
    
    # Полный момент (пружина + демпфер)
    τ = spring.stiffness * Δθ + spring.damping * Δω
    
    # Добавляем в угловые ускорения (делим на инерцию)
    rhs[bd1_vel_dofs[3]] += τ / bd1.inertia
    rhs[bd2_vel_dofs[3]] -= τ / bd2.inertia  # Отрицательный момент на второе тело
end

# Вспомогательная функция для расчета момента пружины
function get_spring_moment(spring::TorsionalSpring, θ1::Float64, θ2::Float64)
    Δθ1 = θ1 - spring.rest_angle
    Δθ2 = θ2 - spring.rest_angle
    return spring.stiffness * (Δθ1 - Δθ2)
end

# Функция для получения энергии пружины
function get_spring_energy(spring::TorsionalSpring, θ1::Float64, θ2::Float64)
    Δθ1 = θ1 - spring.rest_angle
    Δθ2 = θ2 - spring.rest_angle
    Δθ_rel = Δθ1 - Δθ2
    return 0.5 * spring.stiffness * Δθ_rel^2
end