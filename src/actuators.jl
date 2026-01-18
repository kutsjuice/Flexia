mutable struct PositionMotor2D <: AbstractPositionActuator2D
    joint::HingeJoint
    target_angle::Float64
    index::Int64

    function PositionMotor2D(joint::HingeJoint, target_angle::Float64)
        return new(joint, target_angle, -1)
    end
end



function setangle!(act::PositionMotor2D, angle::Float64)
    act.target_angle = angle
    return nothing
end

numberofdofs(::PositionMotor2D) = 1

function get_lms(sys::MBSystem2D, act::PositionMotor2D)
    
    last_lm = sys.lmdofs[act.index]
    return SA[last_lm]
end

function add_to_rhs!(rhs, state, sys::MBSystem2D, act::PositionMotor2D)
    joint = act.joint
    bd1 = joint.body1
    bd2 = joint.body2

    bd1_p_dofs = get_body_position_dofs(sys, bd1)
    bd1_v_dofs = get_body_velocity_dofs(sys, bd1)
    bd2_p_dofs = get_body_position_dofs(sys, bd2)
    bd2_v_dofs = get_body_velocity_dofs(sys, bd2)

    θ1 = state[bd1_p_dofs[3]]
    θ2 = state[bd2_p_dofs[3]]

    ω1 = state[bd1_v_dofs[3]]
    ω2 = state[bd2_v_dofs[3]]

    lms = get_lms(sys, act)
    λ = state[lms[1]]


    # Constraint: θ2 - θ1 - target = 0
    rhs[lms[1]] = θ2 - θ1

    rhs[bd1_p_dofs[3]] -= λ
    rhs[bd2_p_dofs[3]] += λ
end

function propagate_targets!(sys::MBSystem2D, act::PositionMotor2D)
    lms = get_lms(sys, act)
    sys.targets[lms[1]] = act.target_angle
    return nothing
end