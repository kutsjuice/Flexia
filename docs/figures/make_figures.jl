# Generate the figures used in the manual.
#
# Run from the repository root:
#     julia --project=examples docs/figures/make_figures.jl
#
# The output PNG files are committed, so building the documentation does not
# require Julia.

using Flexia
using CairoMakie
using StaticArrays

const FIGDIR = @__DIR__

# ---------------------------------------------------------------------------
# 1. Gallery of body visuals
# ---------------------------------------------------------------------------

function body_visual_gallery()
    labels = ["BarVisual", "RectVisual", "CircleVisual", "Link2DVisual", "MeshVisual"]
    visuals = Any[
        BarVisual(length = 1.0),                       # width = 0 -> a line
        RectVisual(width = 0.9, height = 0.35),
        CircleVisual(radius = 0.3),
        Link2DVisual(length = 0.9, width = 0.28),
        MeshVisual(points = Point2f[(0.0, 0.0), (0.9, 0.0), (1.0, 0.25),
                                    (0.8, 0.5), (0.1, 0.4)]),
    ]
    offsets = [0.0, 1.6, 3.2, 4.8, 6.4]

    fig = Figure(size = (1200, 300), backgroundcolor = :white)
    ax = Axis(fig[1, 1], aspect = DataAspect(),
              limits = (-0.75, 7.15, -0.75, 0.75))
    hidedecorations!(ax)
    hidespines!(ax)

    for (i, v) in enumerate(visuals)
        pts = Flexia.body_contour(offsets[i], 0.0, 0.0, v)
        if v isa BarVisual && v.width <= 0
            lines!(ax, pts; color = :steelblue, linewidth = 4)
        else
            poly!(ax, pts; color = (:steelblue, 0.18), strokecolor = :steelblue,
                  strokewidth = 2)
        end
        text!(ax, offsets[i], -0.5; text = labels[i], align = (:center, :top),
              fontsize = 14)
    end

    save(joinpath(FIGDIR, "body-visuals.png"), fig; px_per_unit = 3)
    return nothing
end

# ---------------------------------------------------------------------------
# 2. Crank-slider mechanism (same structure as examples/crankshaft.jl)
# ---------------------------------------------------------------------------

function crank_slider_figure()
    ground = Body2D(1.0e6, 1.0e6; length = 0.01)
    crank = Body2D(10.0, 1.0; length = 1.0)
    connector = Body2D(5.0, 0.5; length = 2.0)
    slider = Body2D(5.0, 0.1; length = 0.25)

    sys = MBSystem2D()
    add!(sys, ground)
    add!(sys, crank)
    add!(sys, connector)
    add!(sys, slider)

    ground_joint = FixedJoint(ground)
    setposition!(ground_joint, SA[0.0, 0.0])
    setrotation!(ground_joint, 0.0)
    add!(sys, ground_joint)

    crank_hinge = HingeJoint(ground, crank)
    set_position_on_first_body!(crank_hinge, SA[0.0, 0.0])
    set_position_on_second_body!(crank_hinge, SA[-0.5, 0.0])
    add!(sys, crank_hinge)

    rod_hinge = HingeJoint(crank, connector)
    set_position_on_first_body!(rod_hinge, SA[0.5, 0.0])
    set_position_on_second_body!(rod_hinge, SA[-1.0, 0.0])
    add!(sys, rod_hinge)

    slider_hinge = HingeJoint(connector, slider)
    set_position_on_first_body!(slider_hinge, SA[1.0, 0.0])
    set_position_on_second_body!(slider_hinge, SA[0.0, 0.0])
    add!(sys, slider_hinge)

    slider_joint = SliderJoint(ground, slider)
    set_position_on_first_body!(slider_joint, SA[2.3028, 0.0])
    set_position_on_second_body!(slider_joint, SA[0.0, 0.0])
    set_direction_on_first_body!(slider_joint, SA[1.0, 0.0])
    set_direction_on_second_body!(slider_joint, SA[1.0, 0.0])
    add!(sys, slider_joint)

    add!(sys, PositionMotor2D(crank_hinge, 0.0, 0.0))
    add!(sys, LinearSpring(slider_joint; stiffness = 20.0, vis_N = 8))

    assemble!(sys)

    # A bent configuration, consistent with the link lengths:
    # crank angle 60 degrees, connector from the crank end to the slider.
    θ1 = π / 3
    crank_center = SA[0.5 * cos(θ1), 0.5 * sin(θ1), θ1]
    crank_end = SA[cos(θ1), sin(θ1)]
    slider_pos = SA[crank_end[1] + sqrt(4.0 - crank_end[2]^2), 0.0, 0.0]
    connector_center = SA[(crank_end[1] + slider_pos[1]) / 2,
                          crank_end[2] / 2,
                          atan(-crank_end[2], slider_pos[1] - crank_end[1])]

    state = zeros(number_of_dofs(sys))
    for (body, values) in ((crank, crank_center),
                           (connector, connector_center),
                           (slider, slider_pos))
        dofs = get_body_position_dofs(sys, body)
        state[dofs] .= values
    end

    fig, ax, p = systemfigure(sys, reshape(state, :, 1); show_forces = false)
    save(joinpath(FIGDIR, "crank-slider.png"), fig; px_per_unit = 3)
    return nothing
end

body_visual_gallery()
crank_slider_figure()
