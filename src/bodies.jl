mutable struct Body2D <: AbstractBody2D
    mass::Float64
    inertia::Float64

    forces::Vector{Function}
    index::Int64

    # Визуальный компаньон: как тело выглядит. Длина звена переехала сюда
    # (бывшее поле `length`) — в динамике она не участвовала никогда.
    visual::AbstractVisual

    # `length` сохранён как kwarg ради совместимости: он конструирует
    # стандартный BarVisual. Для других форм передавайте `visual = ...`.
    function Body2D(mass, inertia; length = 1.0, visual::AbstractVisual = BarVisual(length = length))
        return new(mass, inertia, [
                (x, t) -> 0,
                (x, t) -> 0,
                (x, t) -> 0
            ], -1, visual)
    end
end
get_index(body::Body2D) = body.index
number_of_dofs(::Body2D) = 3

visual(body::Body2D) = body.visual

function add!(sys::MBSystem2D, body::Body2D)
    push!(sys.bodies, body)
    # sys.bodiesnum += 1
    bodydofs = 2 * number_of_dofs(body)
    lbd = last_body_dof(sys)
    push!(sys.bodiesdofs, lbd + bodydofs)
    body.index = length(sys.bodies)
end

function get_body_position_dofs(sys::MBSystem2D, body::Body2D)
    bid = get_index(body)
    return SA[
        sys.bodiesdofs[bid]-5,
        sys.bodiesdofs[bid]-4,
        sys.bodiesdofs[bid]-3,
    ]

end
function get_body_velocity_dofs(sys::MBSystem2D, body::Body2D)
    bid = get_index(body)
    return SA[
        sys.bodiesdofs[bid]-2,
        sys.bodiesdofs[bid]-1,
        sys.bodiesdofs[bid]-0,
    ]
end

function get_body_generalized_dofs(sys::MBSystem2D, body::Body2D)
    bid = get_index(body)
    lg_dof = Int(sys.bodiesdofs[bid] ÷ 2)
    return SA[
        lg_dof - 2,
        lg_dof - 1,
        lg_dof - 0,
    ]
end

function add_to_rhs!(rhs, state, sys, body)
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

    rhs[velocity_dofs[1]] += body.forces[1](state[[position_dofs; velocity_dofs]], state[end])
    rhs[velocity_dofs[2]] += body.forces[2](state[[position_dofs; velocity_dofs]], state[end])
    rhs[velocity_dofs[3]] += body.forces[3](state[[position_dofs; velocity_dofs]], state[end])
end

"""
    get_boundary_points(sys, body, dofs) -> (Point2f, Point2f)

Совместимость: две крайние точки контура тела. Реальная геометрия — в
[`body_contour`](@ref) (`src/plotdata.jl`), которая диспетчеризуется по
визуальному компаньону тела.
"""
function get_boundary_points(sys::MBSystem2D, body::Body2D, dofs::AbstractVector{Float64})::Tuple{Point2f, Point2f}
    pts = body_contour(sys, body, dofs)
    isempty(pts) && return Point2f(0, 0), Point2f(0, 0)
    return pts[1], pts[end]
end
