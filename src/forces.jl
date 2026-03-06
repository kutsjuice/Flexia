
mutable struct TorsionalSpring <: AbstractForce2D
    body1::Body2D
    body2::Body2D
    stiffness::Float64
    rest_angle::Float64
    damping::Float64
    index::Int64
    
    function TorsionalSpring(body1::Body2D, body2::Body2D, stiffness::Float64=1.0, rest_angle::Float64=0.0, damping::Float64=0.0)
        return new(body1, body2, stiffness, rest_angle, damping, -1)
    end
end

# number_of_dofs(::TorsionalSpring) = 0  # пружина не добавляет лагранжевых множителей

function setid!(force::T, index) where T<:AbstractForce2D
    force.index = index
end

function add!(sys::MBSystem2D, force::AbstractForce2D)
    push!(sys.forces, force)
    setid!(force, length(sys.forces))
end

function add_joint_to_rhs!(rhs, state, sys::MBSystem2D, force::TorsionalSpring)
    bd1 = force.body1
    bd2 = force.body2
    
    bd1_pos_dofs = get_body_position_dofs(sys, bd1)
    bd1_vel_dofs = get_body_velocity_dofs(sys, bd1)
    bd2_pos_dofs = get_body_position_dofs(sys, bd2)
    bd2_vel_dofs = get_body_velocity_dofs(sys, bd2)
    
    θ1 = state[bd1_pos_dofs[3]]
    θ2 = state[bd2_pos_dofs[3]]
    ω1 = state[bd1_vel_dofs[3]]
    ω2 = state[bd2_vel_dofs[3]]
    
    Δθ = θ1 - θ2 - force.rest_angle
    Δω = ω1 - ω2
    
    # Полный момент (пружина + демпфер)
    τ = force.stiffness * Δθ
    + force.stiffness * Δθ^3 
    - force.damping * Δω
    
    # Добавляем в угловые скорости тел
    rhs[bd1_vel_dofs[3]] += -τ 
    rhs[bd2_vel_dofs[3]] += τ  
end

# Вспомогательная функция для расчета момента пружины
function get_spring_moment(force::TorsionalSpring, θ1::Float64, θ2::Float64)
    Δθ1 = θ1 - force.rest_angle
    Δθ2 = θ2 - force.rest_angle
    return force.stiffness * (Δθ1 - Δθ2)
end

# Функция для получения энергии пружины
function get_spring_energy(force::TorsionalSpring, θ1::Float64, θ2::Float64)
    Δθ1 = θ1 - force.rest_angle
    Δθ2 = θ2 - force.rest_angle
    Δθ_rel = Δθ1 - Δθ2
    return 0.5 * force.stiffness * Δθ_rel^2
end

function get_torsionalSpring_point(system::MBSystem2D, force::TorsionalSpring, state::AbstractVector{Float64})
    bd1 = force.body1
    pos_dofs1 = get_body_position_dofs(system, bd1)
    _xi1 = state[pos_dofs1[1]]
    _yi1 = state[pos_dofs1[2]]
    _θi1 = state[pos_dofs1[3]]

    bd2 = force.body2
    pos_dofs2 = get_body_position_dofs(system, bd2)
    _xi2 = state[pos_dofs2[1]]
    _yi2 = state[pos_dofs2[2]]
    _θi2 = state[pos_dofs2[3]]

    _xi = _xi1 + bd1.length*cos(_θi1)
    _yi = _yi1 + bd1.length*sin(_θi1)

    return Point2f(_xi ,_yi)
end