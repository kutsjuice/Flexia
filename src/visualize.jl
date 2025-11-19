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

function Makie.lift(system, solution, joint::TorsionalSpring, i::Observable)
    p =  lift(i) do value
        point = get_torsionalSpring_point(system, joint, view(solution, :, value)) 
        return point;
    end
    return p;
end

function draw!(ax, joint::AbstractJoint2D, system::MBSystem2D, solution, iter::Observable)
end
function draw!(ax, joint::HingeJoint, system::MBSystem2D, solution, iter::Observable)
    hinge_point = lift(system, solution, joint, iter);
    scatter!(ax, hinge_point);
end

function draw!(ax, joint::TorsionalSpring, system::MBSystem2D, solution, iter::Observable)
    hinge_point = lift(system, solution, joint, iter);
    p =  lift(i) do value
        tors_point = get_torsionalSpring_point(system, joint, view(solution, :, value)) 
        return tors_point;
    end
    # Генерируем точки спирали
    θ = 0.:0.05:2*pi
    a = 0.1273
    r = a * θ
    spirla_points = Point2f(tors_point._xi .+ r .* cos.(θ), tors_point._yi .+ r .* sin.(θ))  
    # Рисуем
    lines!(ax, spirla_points, color=:blue, linewidth=2)
    scatter!(ax, hinge_point, color=:red, markersize=8)
end

function animate(sys::MBSystem2D, sol, time_span, filename; framerate=60, limits = (-1, 1, 1, 1))
    fig = Figure()
    iter = Observable(1)
    ax = Axis(fig[1, 1])

    for body in bodies(sys)
        bar = lift(sys, sol, body, iter)
        lines!(ax, bar)
    end

    for joint in joints(sys)
        draw!(ax, joint, sys, sol, iter)
    end

    limits!(ax, limits...)

    record(fig, filename, 1:5:length(time_span);
        framerate=framerate) do t
        iter[] = t
    end
end