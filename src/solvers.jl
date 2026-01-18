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
        sys.prestep(0)
        #cros_step(sys.func, sys.jac, mass, time_step, sol[:, i-1])
        update_targets!(sys)
        u_cur = sol[:, i-1]
        f_i = sys.rhs(u_cur) - get_targets(sys)
        ζ = (mass - (1 + im) / 2 * time_step * sys.jacobian(u_cur)) \ f_i
        sol[:, i] = u_cur + time_step * real.(ζ);
    end
end

function simulate(sys::MBSystem2D, initial::Vector{Float64}, time_span)

    sol = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
    mass = zeros(number_of_dofs(sys), number_of_dofs(sys))
    mass[1:last_body_dof(sys), 1:last_body_dof(sys)] = I(last_body_dof(sys))
    cros!(sol, initial, mass, sys, step(time_span))
    return sol;
end