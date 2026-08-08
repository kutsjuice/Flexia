# Smoke test for src/coordinate_reduction.jl on the parmech mechanism.
# Same geometry as examples/parmech.jl, but with the ground FixedJoint actually
# registered (the example builds it and never calls add!), so that d = 4.

using Pkg; Pkg.activate("examples")
using Flexia
using StaticArrays
using LinearAlgebra
using ForwardDiff

include(joinpath(@__DIR__, "src", "coordinate_reduction.jl"))

const g = 9.81
crank_len = 1.0; conn_len = 1.0; table_heigh = 1.5
H = 1.5; L = 3.0

ground = Body2D(1e6, 1e6; length=0.01)
crank1 = Body2D(1, 1; length=crank_len)
crank2 = Body2D(1, 1; length=crank_len)
crank3 = Body2D(1, 1; length=crank_len)
crank4 = Body2D(1, 1.001; length=crank_len)
connector1 = Body2D(1, 1; length=conn_len)
connector2 = Body2D(1, 1; length=conn_len)
connector3 = Body2D(1, 1; length=conn_len)
connector4 = Body2D(1, 1; length=conn_len)
table1 = Body2D(1, 1; length=table_heigh)
table2 = Body2D(1, 1; length=table_heigh)

for b in (crank1, crank2, crank3, crank4, connector1, connector2, connector3,
          connector4, table1, table2)
    b.forces[2] = (x, t) -> -b.mass * g
end

ground_joint = FixedJoint(ground)
setposition!(ground_joint, SA[0.0, 0.0]); setrotation!(ground_joint, 0.0)

mkhinge(b1, p1, b2, p2) = begin
    j = HingeJoint(b1, b2)
    set_position_on_first_body!(j, p1); set_position_on_second_body!(j, p2); j
end

hinges = [
    mkhinge(ground, SA[-L/2,  H/2], crank1, SA[-crank_len/2, 0.0]),
    mkhinge(ground, SA[-L/2, -H/2], crank2, SA[-crank_len/2, 0.0]),
    mkhinge(ground, SA[ L/2, -H/2], crank3, SA[ crank_len/2, 0.0]),
    mkhinge(ground, SA[ L/2,  H/2], crank4, SA[ crank_len/2, 0.0]),
    mkhinge(crank1, SA[ crank_len/2, 0.0], connector1, SA[-conn_len/2, 0.0]),
    mkhinge(crank2, SA[ crank_len/2, 0.0], connector2, SA[-conn_len/2, 0.0]),
    mkhinge(crank3, SA[-crank_len/2, 0.0], connector3, SA[ conn_len/2, 0.0]),
    mkhinge(crank4, SA[-crank_len/2, 0.0], connector4, SA[ conn_len/2, 0.0]),
    mkhinge(connector1, SA[ conn_len/2, 0.0], table1, SA[ table_heigh/2, 0.0]),
    mkhinge(connector2, SA[ conn_len/2, 0.0], table1, SA[-table_heigh/2, 0.0]),
    mkhinge(connector3, SA[-conn_len/2, 0.0], table2, SA[-table_heigh/2, 0.0]),
    mkhinge(connector4, SA[-conn_len/2, 0.0], table2, SA[ table_heigh/2, 0.0]),
]

slider = SliderJoint(table1, table2)
set_position_on_first_body!(slider, SA[0.0, 0.0])
set_position_on_second_body!(slider, SA[0.0, 0.0])
set_direction_on_first_body!(slider, SA[0.0, -1.0])
set_direction_on_second_body!(slider, SA[0.0, -1.0])

sys = MBSystem2D()
for b in (ground, crank1, crank2, crank3, crank4, connector1, connector2,
          connector3, connector4, table1, table2)
    add!(sys, b)
end
add!(sys, ground_joint)
for j in hinges; add!(sys, j); end
add!(sys, slider)
assemble!(sys)

# ---- analytic starting guess (same as the example) ------------------------
θc = [π/3, π/3, -π/3, -π/3]
p0 = zeros(33)
setp!(body, x, y, θ) = (p0[3body.index-2] = x; p0[3body.index-1] = y; p0[3body.index] = θ)
setp!(ground, 0.0, 0.0, 0.0)
setp!(crank1, -L/2 + crank_len/2*cos(θc[1]),  H/2 + crank_len/2*sin(θc[1]), θc[1])
setp!(crank2, -L/2 + crank_len/2*cos(θc[2]), -H/2 + crank_len/2*sin(θc[2]), θc[2])
setp!(crank3,  L/2 - crank_len/2*cos(θc[3]), -H/2 - crank_len/2*sin(θc[3]), θc[3])
setp!(crank4,  L/2 - crank_len/2*cos(θc[4]),  H/2 - crank_len/2*sin(θc[4]), θc[4])
for (k, (cr, cn, sgn)) in enumerate(((crank1, connector1, 1), (crank2, connector2, 1),
                                     (crank3, connector3, -1), (crank4, connector4, -1)))
    θcn = -θc[k]
    setp!(cn, p0[3cr.index-2] + sgn*(cos(θc[k])*crank_len/2 + cos(θcn)*conn_len/2),
              p0[3cr.index-1] + sgn*(sin(θc[k])*crank_len/2 + sin(θcn)*conn_len/2), θcn)
end
setp!(table1, p0[3connector1.index-2] + cos(-θc[1])*conn_len/2,
              p0[3connector1.index-1] + sin(-θc[1])*conn_len/2 - table_heigh/2, π/2)
setp!(table2, p0[3connector3.index-2] - cos(-θc[3])*conn_len/2,
              p0[3connector3.index-1] - sin(-θc[3])*conn_len/2 + table_heigh/2, π/2)

# ---- reduction ------------------------------------------------------------
cranks = (crank1, crank2, crank3, crank4)
rm = ReducedModel(sys, angle_dofs(cranks))
println("n = $(n_coords(rm)), m = $(n_constraints(rm)), d = $(n_free(rm))")
println("independent = ", rm.independent)

