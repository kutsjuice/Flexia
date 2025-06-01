
mutable struct FixedJoint <: AbstractJoint2D
    body::Body2D
    pos::SVector{2,Float64}
    θ::Float64
    index::Int64

    function FixedJoint(body::Body2D)
        return new(body, [0, 0], 0, -1)
    end
end

function setposition!(joint::FixedJoint, pos)
    joint.pos[1] = pos[1]
    joint.pos[2] = pos[2]
    return nothing
end

function setrotation!(joint::FixedJoint, θ)
    joint.θ = θ
    return nothing
end

function add!(sys::MBSystem2D, joint::FixedJoint)
    push!(sys.joints, joint)
    # sys.jointsnum += 1;
    last_joint_dof = last_lm_dof(sys) + 3
    push!(sys.lmdofs, last_joint_dof)
    joint.index = length(sys.joints)
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

    last_joint_dof = get_lms(sys, joint)

    joint_dofs = [
        last_joint_dof - 2,
        last_joint_dof - 1,
        last_joint_dof,
    ]

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

end

function add_joint_to_rhs!(rhs, state, sys::MBSystem2D, joint::HingeJoint)
    bd1 = joints.body1
    bd2 = joints.body2

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

    rhs[lm[1]] = (_xi + xci * cos(_θi) - yci * sin(_θi)) -
                 (_xj + xcj * cos(_θj) - ycj * sin(_θj))
    rhs[lm[2]] = (_yi + xci * sin(_θi) + yci * cos(_θi)) -
                 (_yj + xcj * sin(_θj) + ycj * sin(_θj))


end