function Makie.lift(system, solution, body::Body2D, i::Observable)
    return lift(i) do value
        points = Vector{Point2f}(undef, 2)
        points[1:2] .= get_boundary_points(system, body, view(solution, :, value))
        return points
    end
end
function Makie.lift(system, solution, joint::AbstractJoint2D, i::Observable)
end

function Makie.lift(system, solution, joint::HingeJoint, i::Observable)
    p =  lift(i) do value
        point = get_hinge_point(system, joint, view(solution, :, value)) 
        return point;
    end
    return p;
end
function draw!(ax, joint::AbstractConnector2D, system::MBSystem2D, solution, iter::Observable)
end
function draw!(ax, joint::HingeJoint, system::MBSystem2D, solution, iter::Observable)
    hinge_point = lift(system, solution, joint, iter);
    scatter!(ax, hinge_point);
end

function animate(sys::MBSystem2D, sol, time_span, filename; framerate=60, limits = (-1, 1, 1, 1))
    fig = Figure()
    iter = Observable(1)
    ax = Axis(fig[1, 1], aspect = DataAspect())

    for body in bodies(sys)
        bar = lift(sys, sol, body, iter)
        lines!(ax, bar)
    end
    for conn in connectors(sys)
        draw!(ax, conn, sys, sol, iter)
    end

    limits!(ax, limits...)

    record(fig, filename, 1:5:length(time_span);
        framerate=framerate) do t
        iter[] = t
    end
end

# ##
# using GLMakie

# fig = Figure()
# ax = Axis(fig[1, 1], aspect = DataAspect())
# t = 0:0.01:2π
# x = sin.(t)
# y = 2cos.(t)
# lines!(ax, x, y)
# valmin = min(minimum(x), minimum(y))
# valmax = max(maximum(x), maximum(y))
# range = valmax - valmin
# xlims!(ax, [valmin - 0.01*range, valmax + 0.01*range])
# ylims!(ax, [valmin - 0.01*range, valmax + 0.01*range])
# fig