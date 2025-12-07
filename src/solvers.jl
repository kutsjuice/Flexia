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

    u = copy(u_cur)
    n = length(u)
    
    history = [copy(u)]
    residuals = [norm(func(u))]

    for iter in max_iter
     Fx = func(u)
        current_residual = norm(Fx)
        push!(residuals, current_residual)
        
        # Проверка сходимости
        if current_residual < tol_e
            return u, iter-1, history, residuals
        end
        
        try
            Ju = jac
            
            # Решаем систему J*dx = -F
            dx = Ju \ (-Fx)
            
            # Обновляем решение
            u .+= dx
            
            push!(history, copy(x))
            
            # Проверка сходимости по изменению решения
            if norm(dx) < tol_e * max(1.0, norm(u))
                return x, iter, history, residuals
            end
            
        catch e
            if isa(e, LinearAlgebra.SingularException)
                @warn "Матрица Якоби вырождена на итерации $iter при x = $x"
                break
            else
                rethrow(e)
            end
        end
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