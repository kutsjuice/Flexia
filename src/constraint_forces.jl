"""
    honest_constraint_forces(sys::MBSystem2D, q, v, v̇) -> Vector{Float64}

Reconstructs the physical constraint/actuator forces (Lagrange multipliers)
at one point of an already-known, consistent trajectory `(q, v, v̇)`, via the
index-1 (acceleration-level) equation of motion

    M v̇ = F(q, v) + Φ_q(q)ᵀ λ

solved here as a linear least-squares problem for λ alone (q, v, v̇ are
supplied externally — e.g. from `simulate`'s own output, or from a planned
spline — not solved for here).

This is deliberately NOT the same λ that `simulate`/`cros!` carries in its
own state: every joint/actuator constraint in this package (`FixedJoint`,
`HingeJoint`, `SliderJoint`, `PositionMotor2D`, `PositionLinearActuator2D`,
...) is enforced at the POSITION level (`Φ(q) = 0`), which makes the overall
system an index-3 DAE with no Baumgarte/GGL stabilization. `q` and `v`
converge properly under `cros!` as the step shrinks, but the λ rows of that
same state do not — they diverge (empirically ~1/h) because the multiplier
has to re-enforce the algebraic constraint from scratch every step. This
function sidesteps that entirely: given a trajectory point that already
satisfies (or nearly satisfies) the constraint, λ is recovered from a
well-posed linear solve, not from integrating the DAE.

Returns one λ per constraint dof, in `sys.lmdofs` order — the same order as
`sol[last_body_dof(sys)+1:last_lm_dof(sys), k]` would be, and indexable with
`get_lms(sys, connector)` exactly as you would that state row.
"""
function honest_constraint_forces(sys::MBSystem2D, q::AbstractVector{<:Real}, v::AbstractVector{<:Real}, v̇::AbstractVector{<:Real})
    @assert sys.assembled "call assemble!(sys) before honest_constraint_forces"

    nb = last_body_dof(sys)
    nλ = last_lm_dof(sys) - nb
    n = number_of_dofs(sys)

    pos_idx = Int[]
    vel_idx = Int[]
    for body in sys.bodies
        append!(pos_idx, get_body_position_dofs(sys, body))
        append!(vel_idx, get_body_velocity_dofs(sys, body))
    end
    @assert length(q) == length(pos_idx)
    @assert length(v) == length(vel_idx)
    @assert length(v̇) == length(vel_idx)

    state = zeros(Float64, n)
    state[pos_idx] = q
    state[vel_idx] = v
    # λ rows (and the trailing time row) are left at zero: every connector's
    # add_to_rhs! reads its own λ out of `state` before adding it into a
    # velocity row (see e.g. actuators.jl), so zeroing them here strips out
    # any constraint-force contribution and leaves the pure applied/coupling
    # force F(q, v) in the velocity rows below.

    Φq = sys.jacobian(state)[(nb+1):(nb+nλ), pos_idx]   # ∂Φ/∂q, n_constraints × n_posdofs
    f = sys.rhs(state)[vel_idx]                          # F(q, v), λ-free

    M = get_mass_matrix(sys)[vel_idx, vel_idx]

    A = Matrix(Φq)'          # n_veldofs × n_constraints (Φ_q(q)ᵀ)
    b = M * collect(v̇) .- f

    return A \ b
end

"""
    honest_constraint_forces(sys::MBSystem2D, sol::AbstractMatrix{<:Real}, dt::Real) -> Matrix{Float64}

Convenience wrapper: applies the point-wise `honest_constraint_forces` at
every column of a `simulate`-style solution `sol` (dofs × time, fixed step
`dt`). Returns an `n_constraints × size(sol, 2)` matrix of reconstructed λ,
in `sys.lmdofs` order.

`v` and `v̇` are estimated by finite differences of `sol`'s POSITION rows
(central, one-sided at the two endpoints) — deliberately NOT read from
`sol`'s own velocity rows. For a body under a rigid `PositionMotor2D`
constraint, `cros!`'s velocity dof stays pinned at ≈0 for the entire
integration regardless of the body's true motion (the same one-stage,
unstabilized index-3 formulation responsible for λ's divergence — see
`honest_constraint_forces(sys, q, v, v̇)` above — apparently also starves the
velocity dof of the implicit correction that would carry true motion into
it), even though the position rows track the true trajectory correctly.
Re-differentiating the (accurate) position history sidesteps that entirely.
"""
function honest_constraint_forces(sys::MBSystem2D, sol::AbstractMatrix{<:Real}, dt::Real)
    nb = last_body_dof(sys)
    nλ = last_lm_dof(sys) - nb
    ntime = size(sol, 2)

    pos_idx = Int[]
    for body in sys.bodies
        append!(pos_idx, get_body_position_dofs(sys, body))
    end
    q = sol[pos_idx, :]   # n_posdofs × ntime

    out = zeros(Float64, nλ, ntime)
    for k in 1:ntime
        v, v̇ = if k == 1
            (q[:, 2] .- q[:, 1]) ./ dt, (q[:, 3] .- 2 .* q[:, 2] .+ q[:, 1]) ./ dt^2
        elseif k == ntime
            (q[:, ntime] .- q[:, ntime-1]) ./ dt, (q[:, ntime] .- 2 .* q[:, ntime-1] .+ q[:, ntime-2]) ./ dt^2
        else
            (q[:, k+1] .- q[:, k-1]) ./ (2dt), (q[:, k+1] .- 2 .* q[:, k] .+ q[:, k-1]) ./ dt^2
        end
        out[:, k] = honest_constraint_forces(sys, q[:, k], v, v̇)
    end
    return out
end
