# Coordinate reduction: index-3 DAE in absolute coordinates -> explicit ODE.
#
#   Flexia:   M_a v̇ = Q(t, p, v) + Cᵀλ,   Φ(p) = 0
#   reduced:  M(q) q̈ + h(q, q̇) = Rᵀ B τ
#
# with  v = R q̇,  v̇ = R q̈ + s,  CR = 0,  M = Rᵀ M_a R,  h = Rᵀ M_a s - Rᵀ Q.
#
# The system must be assembled WITHOUT actuators: in Flexia an actuator is a
# kinematic constraint rather than a source of generalized force, so including
# all of them would leave d = 0 free coordinates. The drive coordinates stay
# free and the torques enter the right-hand side through B.
#
# Index conventions used throughout:
#   * "global" indices address the full state vector of the system
#     (per body: x, y, θ, vx, vy, ω; then the multipliers; then time);
#   * "compact" indices 1:3n_b address the position vector p only, in body
#     order — body k occupies 3k-2, 3k-1, 3k, so 3k is its angle.
# `independent` / `dependent` are always compact indices.

using Flexia
using LinearAlgebra
using ForwardDiff

# ---------------------------------------------------------------------------
# index helpers (compact numbering)
# ---------------------------------------------------------------------------

"Compact indices (x, y, θ) of `body` inside the position vector p."
body_dofs(body::Body2D) = (3body.index - 2, 3body.index - 1, 3body.index)

"Compact index of the angle of `body`."
angle_dof(body::Body2D) = 3body.index

"Compact indices of the angles of `bodies` — the usual choice of independent set."
angle_dofs(bodies) = [angle_dof(b) for b in bodies]

# ---------------------------------------------------------------------------
# model
# ---------------------------------------------------------------------------

struct ReducedModel
    sys::MBSystem2D
    base::Vector{Float64}   # state buffer; only [end] (time) is meaningful here
    pos_cols::Vector{Int}   # global indices of the 3n_b positions
    vel_cols::Vector{Int}   # global indices of the 3n_b velocities
    lm_rows::Vector{Int}    # global indices of the m multipliers
    independent::Vector{Int}
    dependent::Vector{Int}
    Ma::Diagonal{Float64,Vector{Float64}}
end

n_coords(rm::ReducedModel) = length(rm.pos_cols)
n_constraints(rm::ReducedModel) = length(rm.lm_rows)
n_free(rm::ReducedModel) = length(rm.independent)

function _layout(sys::MBSystem2D)
    sys.assembled || error("call assemble!(sys) before building a ReducedModel")
    if any(c -> c isa AbstractActuator2D, sys.connectors)
        error("""the system contains actuators; assemble the model without them.
                 In Flexia an actuator is a kinematic constraint, so its multiplier
                 does not enter the velocity rows and the projection Rᵀ Cᵀ = 0 breaks.""")
    end

    pos_cols = Int[]
    vel_cols = Int[]
    for body in sys.bodies
        append!(pos_cols, get_body_position_dofs(sys, body))
        append!(vel_cols, get_body_velocity_dofs(sys, body))
    end
    lm_rows = collect((last_body_dof(sys)+1):last_lm_dof(sys))
    return pos_cols, vel_cols, lm_rows
end

function ReducedModel(sys::MBSystem2D, independent::AbstractVector{<:Integer}; time::Float64=0.0)
    pos_cols, vel_cols, lm_rows = _layout(sys)

    n = length(pos_cols)
    m = length(lm_rows)
    d = n - m
    d > 0 || error("the system has no free coordinates (3n_b = $n, m = $m)")

    ind = collect(Int, independent)
    allunique(ind) || error("independent coordinates contain duplicates")
    all(i -> 1 <= i <= n, ind) || error("independent coordinates must lie in 1:$n")
    length(ind) == d || error("expected $d independent coordinates (3n_b - m), got $(length(ind))")

    base = zeros(number_of_dofs(sys))
    base[end] = time

    # single source of truth for the inertia: Flexia already assembled it
    Ma = Diagonal([sys.mass[i, i] for i in vel_cols])

    return ReducedModel(sys, base, pos_cols, vel_cols, lm_rows,
                        ind, setdiff(1:n, ind), Ma)
end

"""
    ReducedModel(sys; time=0.0, p0=zeros(3n_b))

Pick the independent set automatically by rank-revealing QR of C(p0): the `m`
columns with the largest pivots become dependent, the rest independent. Use the
explicit constructor when the drive coordinates are known — the automatic choice
is only guaranteed well conditioned near `p0`.
"""
function ReducedModel(sys::MBSystem2D; time::Float64=0.0, p0=nothing)
    pos_cols, _, lm_rows = _layout(sys)
    n = length(pos_cols)
    m = length(lm_rows)
    probe = ReducedModel(sys, collect((m+1):n); time=time)  # placeholder split, indices only
    p = p0 === nothing ? zeros(n) : collect(Float64, p0)
    C = constraint_jacobian(probe, p)
    dependent = sort(qr(C, ColumnNorm()).p[1:m])
    return ReducedModel(sys, setdiff(1:n, dependent); time=time)
end

# ---------------------------------------------------------------------------
# constraints
# ---------------------------------------------------------------------------

"""
    constraint_residual(rm, p)

Φ(p), read straight out of `sys.rhs` (the multiplier rows hold the position-level
constraint residual and depend on positions only). Generic in `eltype(p)`, so it
can be differentiated. Deliberately not `sys.kinematic_constrains!`: that path is
`Vector{Float64}`-only and its `HingeJoint` branch is broken.
"""
function constraint_residual(rm::ReducedModel, p::AbstractVector{T}) where {T}
    st = similar(p, length(rm.base))
    st .= rm.base
    st[rm.pos_cols] .= p
    return rm.sys.rhs(st)[rm.lm_rows]
end

"C = ∂Φ/∂p, size m × 3n_b."
constraint_jacobian(rm::ReducedModel, p::AbstractVector) =
    ForwardDiff.jacobian(pp -> constraint_residual(rm, pp), p)

"""
    assemble_position(rm, q; p0, tol=1e-12, maxiter=50)

Newton solve of Φ(p) = 0 with the independent components pinned to `q`.
The Jacobian is refactorised every iteration (full Newton, not a chord method),
`p0` is the warm start — pass the previous path node.
"""
function assemble_position(rm::ReducedModel, q::AbstractVector{<:Real};
                           p0=zeros(n_coords(rm)), tol::Float64=1e-12, maxiter::Int=50)
    length(q) == n_free(rm) || error("expected $(n_free(rm)) independent coordinates")
    p = collect(Float64, p0)
    p[rm.independent] .= q

    local res
    for _ in 1:maxiter
        Φ = constraint_residual(rm, p)
        res = norm(Φ, Inf)
        res < tol && return p
        F = lu(constraint_jacobian(rm, p)[:, rm.dependent]; check=false)
        (issuccess(F) && pivot_ratio(F) > 1e-12) ||
            error("C_D went singular during assembly (‖Φ‖∞ = $res) — the path " *
                  "crosses a type-II singularity of the chosen independent set")
        p[rm.dependent] .-= F \ Φ
    end
    error("position assembly did not converge in $maxiter iterations (‖Φ‖∞ = $res)")
end

# ---------------------------------------------------------------------------
# reduction
# ---------------------------------------------------------------------------

struct Reduction
    p::Vector{Float64}
    C::Matrix{Float64}
    CD::LU{Float64,Matrix{Float64},Vector{Int}}
    R::Matrix{Float64}
    M::Symmetric{Float64,Matrix{Float64}}
    pivot_ratio::Float64   # see `pivot_ratio`
end

"""
    pivot_ratio(F)

`min|diag U| / max|diag U|` of an LU factorisation — a cheap proxy for how close
`C_D` is to singular, i.e. how close the mechanism is to a type-II singularity
with respect to the chosen independent set. Not a condition number, but it costs
nothing and drops to ~1e-16 exactly when `cond` blows up. Watch it along a path:
a collapse means the independent coordinates stop parameterising the manifold
and the required torques diverge.
"""
function pivot_ratio(F::LU)
    u = abs.(diag(F.U))
    hi = maximum(u)
    return hi == 0 ? 0.0 : minimum(u) / hi
end

"""
    reduce_at(rm, p)

Position-dependent part of the reduction: C, its dependent block factorised,
the null-space basis R (`R_I = I`, `R_D = -C_D⁻¹ C_I`) and `M = Rᵀ M_a R`.
Everything here is reused across the three `h` evaluations needed for the path
coefficients, so build it once per node.
"""
function reduce_at(rm::ReducedModel, p::AbstractVector{Float64}; tol::Float64=1e-12)
    n = n_coords(rm)
    d = n_free(rm)

    C = constraint_jacobian(rm, p)
    F = lu(C[:, rm.dependent]; check=false)
    ratio = issuccess(F) ? pivot_ratio(F) : 0.0
    ratio > tol || error("""C_D is singular at this pose (pivot ratio $(ratio)).
        Either the mechanism is at a type-II singularity with respect to the chosen
        independent coordinates, or the independent set is a bad parameterisation
        here. `ReducedModel(sys; p0=p)` picks a well-conditioned set automatically.""")

    R = zeros(n, d)
    R[rm.independent, :] = I(d)
    R[rm.dependent, :] = -(F \ C[:, rm.independent])

    M = Symmetric(R' * rm.Ma * R)
    return Reduction(collect(Float64, p), C, F, R, M, ratio)
end

"""
    bias_acceleration(rm, red, q̇)

s from `v̇ = R q̈ + s`: `s_I = 0`, `s_D = C_D⁻¹ γ` with `γ = -Ċ v`.

γ is the second directional derivative of Φ along the constraint-consistent
velocity `v = R q̇`, since `d/dδ Φ(p + δv) = C(p + δv) v`. This avoids
differentiating R and avoids assembling Ċ.
"""
function bias_acceleration(rm::ReducedModel, red::Reduction, q̇::AbstractVector{Float64})
    v = red.R * q̇
    γ = -ForwardDiff.derivative(
            ε -> ForwardDiff.derivative(δ -> constraint_residual(rm, red.p .+ δ .* v), ε),
            0.0)
    s = zeros(n_coords(rm))
    s[rm.dependent] = red.CD \ γ
    return s, v
end

"""
    applied_forces(rm, p, v; t)

Q — the generalized applied forces (gravity, springs, body forces) in absolute
coordinates. The multipliers are zeroed so that only the applied part survives.
"""
function applied_forces(rm::ReducedModel, p::AbstractVector{Float64},
                        v::AbstractVector{Float64}; t::Float64=rm.base[end])
    st = copy(rm.base)
    st[rm.pos_cols] .= p
    st[rm.vel_cols] .= v
    st[rm.lm_rows] .= 0.0
    st[end] = t
    return rm.sys.rhs(st)[rm.vel_cols]
end

"h(q, q̇) = Rᵀ M_a s - Rᵀ Q."
function bias_forces(rm::ReducedModel, red::Reduction, q̇::AbstractVector{Float64};
                     t::Float64=rm.base[end])
    s, v = bias_acceleration(rm, red, q̇)
    Q = applied_forces(rm, red.p, v; t=t)
    return red.R' * (rm.Ma * s) - red.R' * Q
end

"""
    reduced_dynamics(rm, p, q̇; t)

`(M, h, R, s, v, Q)` at one state. `M q̈ + h = Rᵀ B τ`.
"""
function reduced_dynamics(rm::ReducedModel, p::AbstractVector{Float64},
                          q̇::AbstractVector{Float64}; t::Float64=rm.base[end])
    red = reduce_at(rm, p)
    s, v = bias_acceleration(rm, red, q̇)
    Q = applied_forces(rm, red.p, v; t=t)
    h = red.R' * (rm.Ma * s) - red.R' * Q
    return (M=red.M, h=h, R=red.R, s=s, v=v, Q=Q, reduction=red)
end

# ---------------------------------------------------------------------------
# actuation
# ---------------------------------------------------------------------------

"""
    actuator_map(rm, pairs)

B, the torque distribution matrix, in compact numbering. Each entry of `pairs`
is `(reaction_body, driven_body)` — bodies, or `nothing` for the reaction body
when the actuator is grounded. Column k carries -1 in the angular row of the
reaction body and +1 in that of the driven body.
"""
function actuator_map(rm::ReducedModel, pairs)
    B = zeros(n_coords(rm), length(pairs))
    for (k, (reaction, driven)) in enumerate(pairs)
        reaction === nothing || (B[angle_dof(reaction), k] = -1.0)
        B[angle_dof(driven), k] = 1.0
    end
    return B
end

"""
    inverse_dynamics(rm, red, q̇, q̈, B; t)

τ from `(Rᵀ B) τ = M q̈ + h`. When every actuator drives an independent angle
against ground, `Rᵀ B = I` and this is just `M q̈ + h`.
"""
function inverse_dynamics(rm::ReducedModel, red::Reduction,
                          q̇::AbstractVector{Float64}, q̈::AbstractVector{Float64},
                          B::AbstractMatrix{Float64}; t::Float64=rm.base[end])
    rhs = red.M * q̈ + bias_forces(rm, red, q̇; t=t)
    return (red.R' * B) \ rhs
end

# ---------------------------------------------------------------------------
# path parameterization
# ---------------------------------------------------------------------------

"""
    path_coefficients(rm, p, q′, q″; t)

τ(s) = a s̈ + b ṡ² + c along a path q(s), with

    a = M q′,   c = h(q, 0),   b = M q″ + h(q, q′) - c.

Valid because h splits into a part quadratic-homogeneous in velocity and a
positional part — which needs Q to be velocity independent (no dampers).
Returns the coefficients in the `M q̈ + h` frame; multiply by `inv(Rᵀ B)` if the
actuation map is not the identity.
"""
function path_coefficients(rm::ReducedModel, p::AbstractVector{Float64},
                           q′::AbstractVector{Float64}, q″::AbstractVector{Float64};
                           t::Float64=rm.base[end])
    red = reduce_at(rm, p)
    c = bias_forces(rm, red, zero(q′); t=t)
    hq = bias_forces(rm, red, q′; t=t)
    a = red.M * q′
    b = red.M * q″ + (hq - c)
    return (a=a, b=b, c=c, M=red.M, R=red.R, reduction=red)
end

# ---------------------------------------------------------------------------
# verification
# ---------------------------------------------------------------------------

"""
    check_reduction(rm, p, q̇; verbose=true)

Sanity checks of the reduction at one state: constraint residual, rank and
conditioning of C, `CR = 0`, positive definiteness of M, the kinetic-energy
identity ½vᵀM_a v = ½q̇ᵀM q̇, and consistency of `s` with `C v̇ = γ`.

Also compares Cᵀ against the multiplier columns of the full Jacobian, row-wise up
to sign. The sign has to be free because Flexia is not consistent about it: for
`HingeJoint`, `FixedJoint` and the second `SliderJoint` row (angle alignment) the
coupling column is `+∂Φ/∂p`, but the first `SliderJoint` row (perpendicular
distance) is distributed as `-∂Φ/∂p` — `joints.jl:335` adds `+λ₁ G_n` to body 1
while `∂Φ₁/∂p₁ = -G_n`. Harmless for the dynamics, since λ just comes out with
the opposite sign and the reduction annihilates the whole term either way; it
only matters if λ is read as a physical reaction force. A naive
`C' == J[vel_cols, lm_rows]` assertion would fail on exactly that row.
"""
function check_reduction(rm::ReducedModel, p::AbstractVector{Float64},
                         q̇::AbstractVector{Float64}; verbose::Bool=true)
    m = n_constraints(rm)
    red = reduce_at(rm, p)
    s, v = bias_acceleration(rm, red, q̇)

    st = copy(rm.base)
    st[rm.pos_cols] .= p
    st[rm.vel_cols] .= v
    J = rm.sys.jacobian(st)
    Ct = J[rm.vel_cols, rm.lm_rows]
    sign_ok = all(1:m) do i
        norm(Ct[:, i] - red.C[i, :]) < 1e-8 || norm(Ct[:, i] + red.C[i, :]) < 1e-8
    end
    flipped = [i for i in 1:m if norm(Ct[:, i] - red.C[i, :]) >= 1e-8]

    γ = -ForwardDiff.derivative(
            ε -> ForwardDiff.derivative(δ -> constraint_residual(rm, p .+ δ .* v), ε), 0.0)

    r = (residual      = norm(constraint_residual(rm, p), Inf),
         rank_C        = rank(red.C),
         m             = m,
         cond_C        = cond(red.C),
         cond_CD       = cond(Matrix(red.C[:, rm.dependent])),
         pivot_ratio   = red.pivot_ratio,
         CR            = norm(red.C * red.R, Inf),
         Cv            = norm(red.C * v, Inf),
         Cs_minus_γ    = norm(red.C * s - γ, Inf),
         posdef_M      = isposdef(red.M),
         energy_error  = abs(0.5 * v' * rm.Ma * v - 0.5 * q̇' * red.M * q̇),
         transpose_ok  = sign_ok,
         flipped_rows  = flipped)

    if verbose
        println("‖Φ‖∞            = ", r.residual)
        println("rank C / m      = ", r.rank_C, " / ", r.m)
        println("cond C          = ", r.cond_C)
        println("cond C_D        = ", r.cond_CD)
        println("pivot ratio     = ", r.pivot_ratio)
        println("‖C R‖∞          = ", r.CR)
        println("‖C v‖∞          = ", r.Cv)
        println("‖C s - γ‖∞      = ", r.Cs_minus_γ)
        println("M posdef        = ", r.posdef_M)
        println("energy mismatch = ", r.energy_error)
        println("Cᵀ match (±)    = ", r.transpose_ok,
                isempty(flipped) ? "" : "  (sign-flipped rows: $flipped)")
    end
    return r
end
