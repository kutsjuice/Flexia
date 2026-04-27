mutable struct MBSystem2D
    bodies::Vector{AbstractBody2D}
    joints::Vector{AbstractJoint2D}
    bodiesdofs::Vector{Int64}
    lmdofs::Vector{Int64}
    assembled::Bool
    has_time::Bool
    time_index::Int64
    rhs::Function
    jacobian::Function
    mass::Matrix{Float64}

    function MBSystem2D()
        default_rhs = (x) -> nothing
        default_jacobian = (x) -> nothing
        return new([], [], [], [], false, false, 0, default_rhs, default_jacobian)
    end
end

function add_time!(sys::MBSystem2D)
    sys.has_time = true
    return nothing
end

number_of_bodies(sys::MBSystem2D) = length(sys.bodies)
bodies(sys::MBSystem2D) = sys.bodies
joints(sys::MBSystem2D) = sys.joints
get_mass_matrix(sys::MBSystem2D) = sys.mass
function last_body_dof(sys::MBSystem2D)
    if (isempty(sys.bodiesdofs))
        return 0
    end
    return (sys.bodiesdofs)[end]
end
function last_lm_dof(sys::MBSystem2D)
    if (isempty(sys.lmdofs))
        return last_body_dof(sys)
    else
        return sys.lmdofs[end]
    end
end

function number_of_dofs(sys) 
    base = last_lm_dof(sys)
    return sys.has_time ? base + 1 : base
end

function assemble!(sys)
    state_length = number_of_dofs(sys)+1

    sys.mass = zeros(state_length, state_length)

    for body in sys.bodies

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

        sys.mass[position_dofs, position_dofs] = I(3)
        sys.mass[velocity_dofs, velocity_dofs] = Diagonal([body.mass, body.mass, body.inertia])
    end 
    sys.mass[end, end] = 1.0

    sys.rhs = (state) -> begin
        ret = similar(state)
        fill!(ret, zero(state[1]))
        for body in sys.bodies
            add_body_to_rhs!(ret, state, sys, body)
            # add_body_to_mass!
        end
        for joint in sys.joints
            add_joint_to_rhs!(ret, state, sys, joint)
        end
        ret[end] = 1.0
        return ret
    end

    sys.jacobian = (state) -> begin
        ret = zeros(state_length, state_length)
        for (i, body) in enumerate(sys.bodies)
            last_dof = sys.bodiesdofs[i]
            position_dofs = [last_dof - 5, last_dof - 4, last_dof - 3]
            velocity_dofs = [last_dof - 2, last_dof - 1, last_dof]
            for p in position_dofs
                ret[p, p] = 1
            end
            for (j, v) in enumerate(velocity_dofs)
                # При наличии времени передаём вектор длины 7
                local_state = sys.has_time ? vcat(state[last_dof-5:last_dof], state[sys.time_index]) : state[last_dof-5:last_dof]
                ret[v, last_dof-5:last_dof] = ForwardDiff.gradient(
                    (x) -> body.forces[j](sys.has_time ? vcat(x, state[sys.time_index]) : x),
                    state[last_dof-5:last_dof])
                ret[v, velocity_dofs] = ForwardDiff.gradient(
                    (x) -> body.forces[j](sys.has_time ? vcat(state[position_dofs], x, state[sys.time_index]) : vcat(state[position_dofs], x)),
                    state[velocity_dofs])
            end
        end
        if sys.has_time
            ret[sys.time_index, :] .= 0.0
        end
        return ret
    end

    sys.assembled = true
end


# function draw!(ax, sys::MBSystem2D, state::Vector)
#     for body in bodies(sys)
#         rod = get_boundary_points(sys, bd1, state)
#         draw!(ax, sys, body, state);
#     end
#     for joint in joints(sys)
#         draw!()
#     end
# end