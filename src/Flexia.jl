# module Flexia
using StaticArrays
# using GLMakie
export MBSystem
export Body2D

mutable struct Body2D
    mass::Float64
    inertia::Float64

    forces::Vector{Function}
    index::Int64
    length::Float64

    function Body2D(mass, inertia; length=1.0)
        return new(mass, inertia, [
                (x) -> 0,
                (x) -> 0,
                (x) -> 0
            ], -1, length)
    end
end

number_of_dofs(::Body2D) = 3

abstract type Joint end
mutable struct FixedJoint <: Joint
    body::Body
    pos::SVector{Float64, 2}
    θ::Float64
    index::Int64

    function FixedJoint(body::Body2D)
        return new(body, [0,0], 0, -1)
    end
end
function setposition!(joint, pos)
    joint.pos[1] = pos[1];
    joint.pos[2] = pos[2];
    return nothing;
end

function setrotation!(joint, θ)
    joint.θ = θ;
    return nothing;
end


mutable struct MBSystem
    bodies::Vector{Body2D}
    # bodiesnum::Int64
    # jointsnum::Int64
    bodiesdofs::Vector{Int64}
    lmdofs::Vector{Int64}
    assembled::Bool
    rhs::Function
    jacobian::Function

    function MBSystem()
        default_rhs = (x) -> nothing
        default_jacobian = (x) -> nothing
        return new([], 0, [], [], false, default_rhs, default_jacobian)
    end
end

function last_body_dof(sys::MBSystem)
    if (isempty(sys.bodiesdofs))
        return 0
    end
    return (sys.bodiesdofs)[end]
end
function last_lm_dof(sys::MBSystem)
    if (isempty(sys.lmdofs))
        return last_body_dof(sys)
    else
        return last_body_dof(sys) + sys.lmdofs[end]
    end
end

number_of_dofs(sys) = last_lm_dof(sys)
function getdofs(sys::MBSystem, body::Body2D)
    if (body.index == -1)
        return 0:0
    elseif (body.index == 1)
        return 1:sys.bodiesdofs[1]
    else
        return (sys.bodiesdofs[body.index-1]+1):sys.bodiesdofs[body.index]
    end
end


function add!(sys::MBSystem, body::Body2D)
    push!(sys.bodies, body)
    # sys.bodiesnum += 1
    bodydofs = 2 * number_of_dofs(body)
    lbd = last_body_dof(sys)
    push!(sys.bodiesdofs, lbd + bodydofs)
    body.index = length(sys.bodies)
end

function add!(sys::MBSystem, joint::FixedJoint)
    push!(sys.joints, joint)
    # sys.jointsnum += 1;
    last_joint_dof = isempty(sys.joints) ? 3 : sys.joints[end] + 3
    push!(sys.lmdofs, last_joint_dof)
    joint.index = length(sys.joints)
end

function add_body_to_rhs!(rhs, state, sys, body)
    last_dof = sys.bodiesdofs[body.index]

    position_dofs = [
        last_dof - 5,
        last_dof - 4,
        last_dof - 3,
    ]

    velocity_dofs = [
        last_dof - 2,
        last_dof - 1,
        last_dof,
    ]
    # assemble velocities
    rhs[position_dofs[1]] = state[velocity_dofs[1]]
    rhs[position_dofs[2]] = state[velocity_dofs[2]]
    rhs[position_dofs[3]] = state[velocity_dofs[3]]

    rhs[velocity_dofs[1]] += body.forces[1](state[position_dofs], state[velocity_dofs]) / body.mass
    rhs[velocity_dofs[2]] += body.forces[2](state[position_dofs], state[velocity_dofs]) / body.mass
    rhs[velocity_dofs[3]] += body.forces[3](state[position_dofs], state[velocity_dofs]) / body.inertia
end

function add_joint_to_rhs!(rhs, state, sys::MBSystem, joint::FixedJoint)
    body = joint.body;
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

    last_joint_dof = last_lm_dof(sys)

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

function assemble!(sys)
    state_length = number_of_dofs(sys)


    sys.rhs = (state) -> begin
        ret = zeros(state_length)
        for  body in sys.bodies
            add_body_to_rhs!(ret, state, sys, body)
        end
        
        return ret
    end

    sys.jacobian = (state) -> begin
        ret = zeros(state_length, state_length)
        for (i, body) in enumerate(sys.bodies)
            last_dof = sys.bodiesdofs[i]

            position_dofs = [
                last_dof - 5,
                last_dof - 4,
                last_dof - 3,
            ]
            for p in position_dofs
                ret[p, p] = 1
            end
            velocity_dofs = [
                last_dof - 2,
                last_dof - 1,
                last_dof,
            ]
            for (j, v) in enumerate(velocity_dofs)
                ret[v, last_dof-5:last_dof] = ForwardDiff().gradient(
                    body.forces[j],
                    state[last_dof-5:last_dof])
                ret[v, velocity_dofs] = ForwardDiff().gradient(
                    (x) -> body.forces[j](state[position_dofs], x),
                    state[velocity_dofs])
            end
        end
        return ret
    end

    assembled = true
end

const g = 9.81

bd1 = Body2D(10, 1)
bd2 = Body2D(10, 1)

bd1.forces[2] = (x) -> -bd1.mass * g
bd2.forces[2] = (x) -> -bd2.mass * g

sys = MBSystem()

add!(sys, bd1)
add!(sys, bd2)

if (!assemble!(sys))
    println("Assembling failed!")
end
# end # module Flexia
func = sys.rhs
using ForwardDiff

# jacoby = (x) -> ForwardDiff.jacobian(func, x)

initial = zeros(12)
func(initial)
# jacoby(initial)