using Flexia
using StaticArrays

const g = 9.81

bd1 = Body2D(10, 1)
bd2 = Body2D(10, 1)
bd3 = Body2D(10, 1)
# bd2 = Body2D(10, 1)

bd1.forces[2] = (x) -> -bd1.mass * g
bd2.forces[2] = (x) -> -bd2.mass * g
bd3.forces[2] = (x) -> -bd2.mass * g


jnt1 = FixedJoint(bd1)
jnt2 = HingeJoint(bd1, bd2)
jnt3 = HingeJoint(bd2, bd3)
set_position_on_second_body!(jnt2, SA[-0.5, 0])

set_position_on_first_body!(jnt3, SA[0.5, 0])
set_position_on_second_body!(jnt3, SA[-0.5, 0])

sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)
add!(sys, bd3)

add!(sys, jnt1)
add!(sys, jnt2)
add!(sys, jnt3)

if (!assemble!(sys))
    println("Assembling failed!")
end
# end # module Flexia
func = sys.rhs



using ForwardDiff

jacoby = (x) -> ForwardDiff.jacobian(func, x)


bd1_x_ind, bd1_y_ind, _ = get_body_position_dofs(sys, bd1)
bd2_x_ind, bd2_y_ind, _ = get_body_position_dofs(sys, bd2)
bd3_x_ind, bd3_y_ind, _ = get_body_position_dofs(sys, bd3)

initial = zeros(number_of_dofs(sys))
initial[bd2_x_ind] = 0.5
func(initial)
jacoby(initial)

mass = zeros(number_of_dofs(sys), number_of_dofs(sys));
for i in 1:last_body_dof(sys)
    mass[i, i] = 1
end
time_span = 0:0.01:10
sol = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
cros!(sol, initial, mass, func, jacoby, step(time_span))
#
using GLMakie

function GLMakie.lift(system, solution, body::Body2D, i::Observable)
    return lift(i) do value
        points = Vector{Point2f}(undef, 2)
        points[1:2] .= get_boundary_points(system, body, solution[:, value])
        return points
    end
end

function animate(sys::MBSystem2D, sol, filename; framerate=60)
    fig = Figure()
    i = Observable(1)
    ax = Axis(fig[1, 1])
    bar1 = lift(sys, sol, bd1, i)
    bar2 = lift(sys, sol, bd2, i)
    bar3 = lift(sys, sol, bd3, i)

    lines!(ax, bar1)
    lines!(ax, bar2)
    lines!(ax, bar3)
    limits!(ax, -4, 4, -4, 4)

    record(fig, filename, 1:5:length(time_span);
        framerate=framerate) do t
        i[] = t
    end
end

animate(sys, sol, "time_animation2.mp4")
