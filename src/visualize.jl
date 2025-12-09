function Makie.lift(system, solution, body::Body2D, i::Observable)
    return lift(i) do value
        points = Vector{Point2f}(undef, 2)
        points[1:2] .= get_boundary_points(system, body, view(solution, :, value))
        return points
    end
end
function Makie.lift(system, solution, joint::AbstractJoint2D, i::Observable)
end

function Makie.lift(system, solution, joint::FixedJoint, i::Observable)
    p =  lift(i) do value
        point = get_fixed_point(system, joint, view(solution, :, value))

        bd = joint.body
        pos_dofs1 = get_body_position_dofs(system, bd)
        _xi1, _yi1, _θi = view(solution, :, value)[pos_dofs1]

        x0 = point[1]
        y0 = point[2]

        R1 = get_lms(system, joint)

        lms11 = 0.001 * solution[R1[1], value]
        lms12 = 0.001 * solution[R1[2], value]
        N1 = 100
        N0 = 4
        N = N1 + N0

        start_angel = _θi
        end_angel = π/2 + π/4

        n1 = 2

        r0 = 0.3 / n1
        r1 = 0.6 / n1
        
        t = LinRange(start_angel, end_angel, N - 1)
        R = LinRange(r0, r1, N - 1)

        points = Vector{Point2f}(undef, N)

        x10 = x0
        y10 = y0
        p10 = Point2f(x10,y10)
        push!(points, p10)

        x11 = x0 + lms11
        y11 = y0 + lms12
        p11 = Point2f(x11,y11)
        push!(points, p11)
        
        x20 = x0
        y20 = y0
        p20 = Point2f(x20, y20)
        push!(points, p20)

        for j in 1:N1
            x3 = R[j] * cos(t[j]) + x0  
            y3 = R[j] * sin(t[j]) + y0
            p3 = Point2f(x3,y3)
            push!(points, p3)
        end
        
        p_end = Point2f(x0, y0)
        push!(points, p_end)

        return points;
    end
    return p;
end

function Makie.lift(system, solution, joint::HingeJoint, i::Observable)
    p =  lift(i) do value
        point = get_hinge_point(system, joint, view(solution, :, value)) 
        return point;
    end
    return p;
end

# function Makie.lift(system, solution, joint::TorsionalSpring, i::Observable)
#     p =  lift(i) do value
#         a = 0.1273
#         r = a * value
#         point = get_torsionalSpring_point(system, joint, view(solution, :, value))
#         point2 = Point2f(point[1] .+ r .* cos.(value), point[2] + .+ r .* sin.(value) )
#         return point2;
#     end
#     return p;
# end
function Makie.lift(system, solution, joint::TorsionalSpring, i::Observable)
    p = lift(i) do value
        point = get_torsionalSpring_point(system, joint, view(solution, :, value)) 
        bd1 = joint.body1
        pos_dofs1 = get_body_position_dofs(system, bd1)
        _xi1, _yi1, _θi1 = view(solution, :, value)[pos_dofs1]

        bd2 = joint.body2
        pos_dofs2 = get_body_position_dofs(system, bd2)
        _xi2, _yi2, _θi2 = view(solution, :, value)[pos_dofs2]

        start_angel = _θi1 + π
        end_angel = _θi2 + 4*π

        n1 = 2

        r0 = 0.3 / n1
        r1 = 0.6 / n1
        N = 100

        t = LinRange(start_angel, end_angel, N)
        R = LinRange(r0, r1, N)
        x0 = point[1]
        y0 = point[2]

        points = Vector{Point2f}(undef, N)
        for j in 1:N
            x = R[j] * cos(t[j]) + x0  
            y = R[j] * sin(t[j]) + y0
            points[j] = Point2f(x,y)
        end
        return points;
    end
    return p
end

function draw!(ax, joint::AbstractJoint2D, system::MBSystem2D, solution, iter::Observable)
end
function draw!(ax, joint::HingeJoint, system::MBSystem2D, solution, iter::Observable)
    hinge_point = lift(system, solution, joint, iter);
    scatter!(ax, hinge_point);
end

function draw!(ax, joint::TorsionalSpring, system::MBSystem2D, solution, iter::Observable)
    hinge_point = lift(system, solution, joint, iter);
    lines!(ax, hinge_point);
    # stati = UndefInitializer{Float64}[undef,3]
    # bd1 = joint.body1
    # pos_dofs1 = get_body_position_dofs(system, bd1)
    # _xi1, _yi1, _θi1 = stati[pos_dofs1]

    # bd2 = joint.body2
    # pos_dofs2 = get_body_position_dofs(system, bd2)
    # _xi2, _yi2, _θi2 = stati[pos_dofs2]

    # _delta_θ = _θi1:1:_θi2

#     circle_with_hole = BezierPath([
#     MoveTo(Point(1, 0)),
#     EllipticalArc(Point(0, 0), 1, 1, 0, 0, 2pi),
#     MoveTo(Point(0.5, 0.5)),
#     LineTo(Point(0.5, -0.5)),
#     LineTo(Point(-0.5, -0.5)),
#     LineTo(Point(-0.5, 0.5)),
#     ClosePath(),
# ])
#     spiral = BezierPath([
#     MoveTo(Point(0, 0)),
#     LineTo(Point(5.07, 3.197)),
#     LineTo(Point(2.54, 7.05)),
#     LineTo(Point(-9.51, 5.12)),
#     LineTo(Point(-12.11, -5.567)),
#     LineTo(Point(0., -8.841)),
#     LineTo(Point(13.32, 0.)),
#     LineTo(Point(5.077, 13.438)),
#     LineTo(Point(-9.509, 13.438)),
#     LineTo(Point(-9.509, 15.438)),
#     LineTo(Point(5.489, 15.273)),
#     LineTo(Point(14.991, 0.)),
#     LineTo(Point(0., -10.645)),
#     LineTo(Point(-13.33, -6.376)),
#     LineTo(Point(-10.811, 6.08)),
#     LineTo(Point(2.538, 8.558)),
#     LineTo(Point(6.658, 3.197)),
#     LineTo(Point(0.593, -2.019)),
#     ClosePath(),
# ])

    # scatter!(ax, hinge_point, marker = spiral, color=:red,rotation = range(_θi1, _θi2, length = 6)[1:end-1], markersize=1);

    # spiral = lift(system, solution, joint, iter)
    # lines!(ax, spiral)
end

function draw!(ax, joint::FixedJoint, system::MBSystem2D, solution, iter::Observable)
    fixed_point = lift(system, solution, joint, iter);
    lines!(ax, fixed_point);
end

function animate(sys::MBSystem2D, sol, time_span, filename; framerate=60, limits = (-1, 1, 1, 1))
    fig = Figure()
    iter = Observable(1)
    ax = Axis(fig[1, 1], aspect = DataAspect())

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