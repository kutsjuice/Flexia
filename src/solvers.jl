function cros_step(func::Function, jac::Function, mass::Matrix{T}, time_step::T, u_cur::Vector{T})::Vector{T} where T<:Real
    ζ = (mass - (1 + im) / 2 * time_step * jac(u_cur)) \ func(u_cur)

    return u_cur + time_step * real.(ζ)
end



function cros!(sol::Matrix{T}, u0::Vector{T}, mass::Matrix{T}, func::Function, jac::Function, time_step::T) where T<:Real
    @assert size(sol, 1) == size(u0, 1)
    @assert size(sol, 1) == size(func(u0), 1)
    sol[:, 1] = u0
    for i in 2:size(sol, 2)
        sol[:, i] = cros_step(func, jac, mass, time_step, sol[:, i-1])
    end
end


function cros!(sol::Matrix{T}, u0::Vector{T}, mass::Matrix{T}, sys::MBSystem2D, time_step::T) where T<:Real
    @assert size(sol, 1) == size(u0, 1)
    @assert size(sol, 1) == size(sys.rhs(u0), 1)
    sol[:, 1] = u0
    for i in 2:size(sol, 2)
        u_cur = sol[:, i-1]
        sys.prestep(u_cur)
        update_targets!(sys)
        f_i = sys.rhs(u_cur)
        ζ = (mass - (1 + im) / 2 * time_step * sys.jacobian(u_cur)) \ f_i
        sol[:, i] = u_cur + time_step * real.(ζ);
    end
end

function simulate(sys::MBSystem2D, initial::Vector{Float64}, time_span)

    sol = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
    mass = sys.mass
    cros!(sol, initial, mass, sys, step(time_span))
    return sol;
end