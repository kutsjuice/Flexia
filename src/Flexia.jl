# module Flexia
using StaticArrays
# using GLMakie
export MBSystem
export Body2D

struct Body2D
    mass::Float64
    inertia::Float64

    forces::Vector{Function}

    function Body2D(mass, inertia)
        return new(mass, inertia, [
            (x) -> 0,
            (x) -> 0,
            (x) -> 0
        ])
    end
end

number_of_dofs(::Body2D) = 3

mutable struct MBSystem
    bodies::Vector{Body2D}
    bodiesnum::Int64
    bodiesdofs::Vector{Int64}
    assebmled::Bool
    rhs::Function
    jacobian::Function

    function MBSystem()
        default_rhs = (x) -> nothing

        return new([], 0, [], false, default_rhs)
    end
end

function last_body_dof(sys::MBSystem)
    if (isempty(sys.bodiesdofs))
        return 0
    end
    return (sys.bodiesdofs)[end]
end

number_of_dofs(sys) = last_body_dof(sys)

function add!(sys::MBSystem, body::Body2D)
    push!(sys.bodies, body)
    sys.bodiesnum += 1
    bodydofs = 2 * number_of_dofs(body)
    lbd = last_body_dof(sys)
    push!(sys.bodiesdofs, lbd + bodydofs)
end




function assemble!(sys)
    state_length = number_of_dofs(sys)


    sys.rhs = (state) -> begin
        ret = zeros(state_length)
        for (i, body) in enumerate(sys.bodies)
            last_dof = sys.bodiesdofs[i]

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
            ret[position_dofs[1]] = state[velocity_dofs[1]]
            ret[position_dofs[2]] = state[velocity_dofs[2]]
            ret[position_dofs[3]] = state[velocity_dofs[3]]

            ret[velocity_dofs[1]] += body.forces[1](state[position_dofs], state[velocity_dofs]) / body.mass
            ret[velocity_dofs[2]] += body.forces[2](state[position_dofs], state[velocity_dofs]) / body.mass
            ret[velocity_dofs[3]] += body.forces[3](state[position_dofs], state[velocity_dofs]) / body.inertia
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
            for (j,v) in enumerate(velocity_dofs)
                ret[v, last_dof-5:last_dof] = ForwardDiff().gradient(
                    body.forces[j], 
                    state[last_dof-5:last_dof]);
                ret[v, velocity_dofs] = ForwardDiff().gradient(
                    (x) -> body.forces[j](state[position_dofs], x), 
                    state[velocity_dofs]);
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