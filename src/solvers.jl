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

function cros!(sol::Matrix{T}, u0::Vector{T}, mass::Matrix{T}, sys::MBSystem2D, time_step::T) where T<:Real
    @assert size(sol, 1) == size(u0, 1)
    @assert size(sol, 1) == size(sys.rhs(u0), 1)
    act_jac = zeros(size(sys.rhs(u0), 1), size(sys.rhs(u0), 1))
    sol[:, 1] = u0
    for i in 2:size(sol, 2)
        #cros_step(sys.func, sys.jac, mass, time_step, sol[:, i-1])
        u_cur = sol[:, i-1]
        sys.prestep(u_cur)
        update_targets!(sys)
        f_i = sys.rhs(u_cur)
        update_act_jacobian!(sys, act_jac)
        jac_ful = sys.jacobian(u_cur) + act_jac
        ζ = (mass - (1 + im) / 2 * time_step * jac_ful) \ f_i
        sol[:, i] = u_cur + time_step * real.(ζ);
    end
end

function simulate(sys::MBSystem2D, initial::Vector{Float64}, time_span)

    sol = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
    mass = sys.mass
    cros!(sol, initial, mass, sys, step(time_span))
    return sol;
end

function vector_norm(x::AbstractVector, p::Real=2)
    if p == 1
        return sum(abs, x)
    elseif p == 2
        return sqrt(sum(abs2, x))
    elseif p == Inf
        return maximum(abs, x)
    elseif p == -Inf
        return minimum(abs, x)
    elseif p > 0
        return (sum(abs(x).^p))^(1/p)
    else
        throw(ArgumentError("Параметр p должен быть положительным или Inf"))
    end
end

function matrix_norm(A::AbstractMatrix, norm_type::Symbol=:frobenius)
    if norm_type === :one
        # 1-норма матрицы (максимальная сумма по столбцам)
        return maximum(sum(abs, A, dims=1)[:])
    elseif norm_type === :inf
        # ∞-норма матрицы (максимальная сумма по строкам)
        return maximum(sum(abs, A, dims=2)[:])
    elseif norm_type === :frobenius
        # Фробениусова норма
        return sqrt(sum(abs2, A))
    elseif norm_type === :two
        # 2-норма матрицы (спектральная норма)
        return svd(A).S[1]  # Наибольшее сингулярное значение
    elseif norm_type === :max
        # Максимальная абсолютная величина элемента
        return maximum(abs, A)
    elseif norm_type === :nuclear
        # Ядерная норма (сумма сингулярных значений)
        return sum(svd(A).S)
    else
        throw(ArgumentError("Неизвестный тип нормы: $norm_type. Допустимые значения: :one, :two, :inf, :frobenius, :max, :nuclear"))
    end
end

function norm(A::AbstractVector, p::Real=2)
    return vector_norm(A, p)
end

function norm(A::AbstractMatrix, norm_type::Symbol=:frobenius)
    return matrix_norm(A, norm_type)
end

# ===========================================================================
# Статика: решение F(u) = 0 методом Ньютона
# ===========================================================================
#
# ПРО СЛОТ ВРЕМЕНИ. `number_of_dofs(sys) = last_lm_dof + 1`: последняя
# компонента состояния — время, и `sys.rhs` возвращает в этой строке константу
# `ret[end] = 1.0`. В статике время не неизвестное, а ПАРАМЕТР, поэтому оно
# исключается из решаемой системы: `newton_step` отрезает последнюю компоненту,
# подставляет её значение как параметр и решает систему размерности n−1.
#
# Раньше решалась полная система n×n, а её якобиан структурно вырожден
# (нулевая строка времени) — `J \ (-F)` бросал `SingularException`, исключение
# глушил пустой `catch`, и решатель молча возвращал начальное приближение.
#
# Уравнения статики (абсолютные координаты, Лагранж I рода):
#     Q(q) − Φ_qᵀ λ = 0,    Φ(q) = 0,
# а строки `rhs[vel] = ...` при обращении в ноль сами дают v = 0.

"""
    newton_solve(F, J, u0; max_iter, tol, tol_step, max_backtrack)
        -> (u, status, iterations, residuals)

Метод Ньютона для `F(u) = 0` в чистом виде (без слота времени): `u0` — полный
вектор неизвестных, `F` и `J` — невязка и её якобиан той же размерности.

Статус:

  * `:converged` — `‖F‖ < tol`;
  * `:stalled` — шаг стал меньше `tol_step`, а `‖F‖` так и не упала ниже `tol`;
  * `:max_iter` — исчерпан `max_iter`.

Шаг демпфируется по невязке: полный шаг Ньютона принимается, только если он
уменьшает `‖F‖`, иначе дробится пополам (до `max_backtrack` раз). Исключения НЕ
глушатся: вырожденный якобиан пробрасывается наружу как `SingularException`.
"""
function newton_solve(F::Function, J::Function, u0::AbstractVector{T};
                      max_iter::Integer = 100,
                      tol::Real = 1e-10,
                      tol_step::Real = 1e-12,
                      max_backtrack::Integer = 20) where T <: Real

    u = collect(T, u0)
    residuals = Float64[norm(F(u))]

    residuals[1] < tol && return u, :converged, 0, residuals

    for iter in 1:max_iter
        Fu = F(u)
        r = norm(Fu)

        du = -(J(u) \ Fu)

        # Демпфирование: первый alpha, который уменьшает невязку.
        alpha = one(T)
        u_new = u + alpha .* du
        for _ in 1:max_backtrack
            norm(F(u_new)) < r && break
            alpha /= 2
            u_new = u + alpha .* du
        end

        u = u_new
        push!(residuals, norm(F(u)))

        residuals[end] < tol && return u, :converged, iter, residuals
        alpha * norm(du) < tol_step && return u, :stalled, iter, residuals
    end

    return u, :max_iter, max_iter, residuals
end

"""
    newton_step(func, jac, u_cur; kwargs...) -> (u, status, iterations, residuals)

Статический шаг Ньютона для системы Flexia: `func`/`jac` — `sys.rhs` и её
якобиан на полном векторе состояния, последняя компонента которого — время.
Время фиксируется на значении `u_cur[end]` и в число неизвестных не входит
(см. комментарий выше); фактическое решение — [`newton_solve`](@ref).

Возвращает полный вектор состояния (со слотом времени на месте `u_cur[end]`)
и статус сходимости.
"""
function newton_step(func::Function, jac::Function, u_cur::AbstractVector{T};
                     max_iter::Integer = 100,
                     tol::Real = 1e-10,
                     tol_step::Real = 1e-12,
                     max_backtrack::Integer = 20) where T <: Real

    full = collect(T, u_cur)
    n = length(full)
    @assert n >= 2 "в состоянии должны быть неизвестные и слот времени"
    t0 = full[end]
    m = n - 1

    F = (u) -> func(vcat(u, t0))[1:m]
    J = (u) -> jac(vcat(u, t0))[1:m, 1:m]

    u, status, iters, residuals = newton_solve(F, J, @view(full[1:m]);
        max_iter = max_iter, tol = tol, tol_step = tol_step,
        max_backtrack = max_backtrack)

    return vcat(u, t0), status, iters, residuals
end

"""
    static_solver!(sol, u0, func, jac; kwargs...) -> NamedTuple

Статический расчёт: находит равновесие `F(u) = 0` из начального приближения
`u0` и записывает его **во все столбцы** `sol`.

Матрица из многих столбцов оставлена для совместимости со старыми вызовами
(`sol[:, end]` — ответ). Раньше число столбцов означало число шагов
релаксации, и точность зависела от него; теперь статическое решение одно и
ищется до сходимости, а число столбцов на результат не влияет.

Возвращает `(status, iterations, residual)`; `status = :converged` означает,
что невязка дошла до `tol`.
"""
function static_solver!(sol::AbstractMatrix{T}, u0::AbstractVector{T},
                        func::Function, jac::Function;
                        max_iter::Integer = 100,
                        tol::Real = 1e-10,
                        tol_step::Real = 1e-12,
                        max_backtrack::Integer = 20) where T <: Real

    @assert size(sol, 1) == length(u0) "размерность sol и u0 не совпадают"
    @assert size(sol, 1) == length(func(u0)) "размерность sol и RHS не совпадают"
    @assert size(sol, 2) >= 1 "нужен хотя бы один столбец"

    u, status, iters, residuals = newton_step(func, jac, u0;
        max_iter = max_iter, tol = tol, tol_step = tol_step,
        max_backtrack = max_backtrack)

    for i in 1:size(sol, 2)
        sol[:, i] = u
    end

    if status !== :converged
        @warn "static_solver!: статический расчёт не сошёлся" status = status iterations = iters residual = residuals[end]
    end

    return (status = status, iterations = iters, residual = residuals[end])
end

#

function find_natural_freqs(sys::MBSystem2D, position)

    bodies_in_order = sys.bodies            # берём из sys, чтобы точно знать порядок
    nb   = last_body_dof(sys)               # = 6 * число тел
    nλ   = last_lm_dof(sys) - nb
    ntot = nb + nλ

    pos_idx = Int[]
    vel_idx = Int[]
    for bd in bodies_in_order
        ix, iy, it = get_body_position_dofs(sys, bd)
        append!(pos_idx, (ix, iy, it))
        vx, vy, wt = get_body_velocity_dofs(sys, bd)
        append!(vel_idx, (vx, vy, wt))
    end
    λ_idx = (nb+1):ntot

    npos = length(pos_idx)

    func = sys.rhs

    jacoby = (x) -> ForwardDiff.jacobian(func, x)


    J = jacoby(position)

    E = Flexia.get_mass_matrix(sys)

    # Блоки якобиана силовой части
    mK = -J[vel_idx, pos_idx]   # = -K_eff / m (включая геом. члены от λ)
    mM = E[vel_idx, vel_idx]

    # Матрица связей (однозначно)
    C = J[λ_idx, pos_idx]

    F_qr = qr(Matrix(C'))

    Qfull = F_qr.Q * Matrix(I, npos, npos)    # принудительно расширяем до 24×24
    r = rank(C)
    V_n = Qfull[:, r+1:end]


    M̃ = V_n' * mM * V_n
    K̃ = V_n' * mK    * V_n

    # Симметризуем
    M̃_sym = Symmetric((M̃ + M̃') / 2)
    K̃_sym = Symmetric((K̃ + K̃') / 2)

    eig = eigen(Matrix(K̃_sym), Matrix(M̃_sym))
    ω² = real.(eig.values)
    order = sortperm(ω²)
    ω²    = ω²[order]
    return sqrt.(ω²) / (2π)
end