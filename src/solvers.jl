function cros_step(func::Function,jac::Function, mass::Matrix{T}, time_step::T, u_cur::Vector{T})::Vector{T} where T <: Real
    ζ = (mass - (1+im)/2 * time_step * jac(u_cur))\func(u_cur);

    return u_cur + time_step * real.(ζ)
end

function cros!(sol::Matrix{T}, u0::Vector{T}, mass::Matrix{T}, func::Function, jac::Function, time_step::T) where T<: Real
    @assert size(sol, 1) == size(u0, 1)
    @assert size(sol, 1) == size(func(u0),1)
    sol[:, 1] = u0
    for i in 2:size(sol, 2)
        sol[:, i] = cros_step(func, jac, mass, time_step, sol[:, i-1])
    end
end

function newton_step(func::Function, jac::Function, u_cur::Vector{T}, max_iter = 1000, tol_e = 1e-5) where T<: Real

    u = u_cur
    history = [u]

    for i in max_iter
         f = func(u)
         f_diff = jac(u)

         if abs(f_diff) < tol_e
            @warn "Метод может не сойтись"
            break
         end
        u_new = u - f / f_diff

        push!(history, u_new)

        if abs(u_new - u) < tol_e || abs(func(u_new)) < tol_e
            return u_new, i, history 
        end

        u = u_new
    end
    
    @warn "Метод не сошёлся"

    return u, max_iter, history
end

function static_solver!(sol::Matrix{T}, u0::Vector{T}, func::Function, jac::Function) where T<: Real
    @assert size(sol, 1) == size(u0, 1)
    @assert size(sol, 1) == size(func(u0),1)
    sol[:, 1] = u0
    
    for i in 2:size(sol, 2)
        sol[:, i] = newton_step(func, jac, sol[:, i-1])
    end

end