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
# Статика: равновесие F(u) = 0
# ===========================================================================
#
# ПРО СЛОТ ВРЕМЕНИ. `number_of_dofs(sys) = last_lm_dof + 1`: последняя
# компонента состояния — время, и `sys.rhs` возвращает в этой строке константу
# `ret[end] = 1.0`. В статике время не неизвестное, а ПАРАМЕТР: оно вырезается
# из решаемой системы (размерность n−1) и берётся либо из `u0[end]`, либо из
# списка `times` — это псевдовремя нагрузки, аналог Load Step в ANSYS Static
# Structural.
#
# Полную систему n×n решать нельзя: её якобиан структурно вырожден (нулевая
# строка времени), `J \ (-F)` бросает `SingularException`. Раньше это
# исключение глушил пустой `catch`, и решатель молча возвращал начальное
# приближение — то есть был no-op.
#
# ДВА РЕЖИМА (см. `static_solver!`): один подшаг (по умолчанию, на выходе
# матрица-столбец) и путь по псевдовремени (несколько столбцов = наборы
# результатов, годятся для анимации нагружения).
#
# СОЛВЕР ПОДКЛЮЧАЕМЫЙ. Контракт:
#     solver(F, J, u0; solver_options...) -> (u, status, iterations, residuals)
# Стандартный — `_newton_solve` (без внешних зависимостей). Другой, например
# из NLsolve.jl, подключается ключом `solver`: библиотека от него не зависит.
#
# Уравнения статики (абсолютные координаты, Лагранж I рода):
#     Q(q) − Φ_qᵀ λ = 0,    Φ(q) = 0,
# а строки `rhs[vel] = ...` при обращении в ноль сами дают v = 0.

"""
    _newton_solve(F, J, u0; max_iter, tol, tol_step, max_backtrack)
        -> (u, status, iterations, residuals)

Стандартный солвер статики Flexia: метод Ньютона для `F(u) = 0` в чистом виде
(без слота времени): `u0` — полный вектор неизвестных, `F` и `J` — невязка и
её якобиан той же размерности.

Статус:

  * `:converged` — `‖F‖ < tol`;
  * `:stalled` — шаг стал меньше `tol_step`, а `‖F‖` так и не упала ниже `tol`;
  * `:max_iter` — исчерпан `max_iter`.

Шаг демпфируется по невязке: полный шаг Ньютона принимается, только если он
уменьшает `‖F‖`, иначе дробится пополам (до `max_backtrack` раз). Исключения НЕ
глушатся: вырожденный якобиан пробрасывается наружу как `SingularException`.

Внешний солвер (например, из NLsolve.jl) должен иметь ту же сигнатуру и
возвращать ту же четвёрку — он подключается ключом `solver` в
[`static_solver!`](@ref).
"""
function _newton_solve(F::Function, J::Function, u0::AbstractVector{T};
                      max_iter::Integer = 100,
                      tol::Real = 1e-10,
                      tol_step::Real = 1e-12,
                      max_backtrack::Integer = 20) where T <: Real

    u = collect(T, u0)
    residuals = Float64[norm(F(u))]

    residuals[1] < tol && return u, :converged, 0, residuals

    for iter in 1:max_iter
        Fu = F(u)
    
        du = -(J(u) \ Fu)

        # Демпфирование: первый alpha, который уменьшает невязку.
        u_new = u + du
        u = u_new
        push!(residuals, norm(F(u)))
        residuals[end] < tol && return u, :converged, iter, residuals
        norm(du) < tol_step && return u, :stalled, iter, residuals
    end

    return u, :max_iter, max_iter, residuals
end

"""
    static_solver!(sol, u0, func, jac; times, solver, solver_options,
                   cutback, min_step, max_cutbacks) -> Vector{NamedTuple}

Статический расчёт: равновесие `F(u) = 0`, где `func`/`jac` — `sys.rhs` и её
якобиан на полном векторе состояния (последняя компонента — время).

Время — не неизвестное, а параметр: оно вырезается из решаемой системы и
берётся из `u0[end]` либо из списка `times` (см. комментарий выше).

Режимы:

  * `times = nothing` (по умолчанию) — ОДНА решённая подзадача при
    `t = u0[end]`: на выходе матрица-столбец. Точность определяется
    сходимостью солвера, а не числом столбцов `sol`.
  * `times = [t₁, t₂, …]` — путь по псевдовремени (аналог Load Step /
    Substeps в ANSYS Static Structural): на каждом `t_k` решается своё
    равновесие, тёплый старт от предыдущего. Столбец `k` — набор результатов
    при `t_k` (в строке `end` лежит само `t_k`), поэтому матрицу можно
    анимировать. Именно так задаётся ramp нагрузки: `force = t -> F_max*t/T`.

`sol`: первые `length(times)` столбцов — результаты подшагов; если столбцов
больше, хвост заполняется последним кадром (`sol[:, end]` — ответ); если
столбцов меньше, чем времён, — `ArgumentError`.

Подключаемый солвер: `solver(F, J, u0; solver_options...) -> (u, status,
iterations, residuals)`; по умолчанию [`_newton_solve`](@ref). Параметры
солвера передаются в `solver_options`, например
`(max_iter = 50, tol = 1e-12)`.

`cutback < 1` включает автоматическое дробление подшага (ANSYS «automatic time
stepping»): если к очередному `t` солвер не сошёлся, интервал делится — в
очередь вставляется промежуточный момент `t + dt·cutback`, и он тоже
становится кадром (кадров может оказаться больше, чем `length(times)`; за
размером `sol` следит вызывающий). Дробление прекращается, когда шаг уходит
ниже `min_step` или когда неудач подряд больше `max_cutbacks` — тогда `@warn`
и остановка.

Возвращает журнал по кадрам:
`(frame, time, status, iterations, residual, step)`.
"""
function static_solver!(sol::AbstractMatrix{T}, u0::AbstractVector{T},
                        func::Function, jac::Function;
                        times = nothing,
                        solver::Function = _newton_solve,
                        solver_options::NamedTuple = NamedTuple(),
                        cutback::Real = 1.0,
                        min_step::Real = 1e-4,
                        max_cutbacks::Integer = 20) where T <: Real

    n = length(u0)
    @assert n >= 2 "в состоянии должны быть неизвестные и слот времени"
    @assert size(sol, 1) == n "размерность sol и u0 не совпадают"
    @assert size(sol, 1) == length(func(u0)) "размерность sol и RHS не совпадают"
    @assert 0 < cutback <= 1 "cutback должен лежать в (0, 1]"

    m = n - 1
    queue = Float64.(times === nothing ? [u0[end]] : times)
    @assert !isempty(queue) "список времён пуст"
    @assert length(queue) == 1 || all(diff(queue) .> 0) "times должны возрастать"
    @assert all(tt -> tt >= u0[end] - 10 * eps(Float64), queue) "times должны быть не раньше u0[end]: псевдовремя идёт только вперёд"
    @assert size(sol, 2) >= length(queue) "в sol меньше столбцов, чем времён в times"

    u = collect(T, u0)
    log = NamedTuple[]
    frame = 1
    sol[:, frame] = u0
    t = u0[end]
    cutbacks = 0

    while !isempty(queue)
        t_target = first(queue)
        if t_target < t - 10 * eps(Float64)
            popfirst!(queue)          # цель уже пройдена после дробления
            continue
        end
        dt = t_target - t

        # Время — параметр: подставляем t_target и решаем систему n−1.
        F = (x) -> func(vcat(x, t_target))[1:m]
        J = (x) -> jac(vcat(x, t_target))[1:m, 1:m]
        u_new, status, iters, residuals =
            solver(F, J, copy(sol[1:m, frame]); solver_options...)

        if status === :converged
            if frame + 1 > size(sol, 2)
                @warn "static_solver!: кадров больше, чем столбцов sol" frames = frame + 1 time = t_target
                break
            end
            frame += 1
            u[1:m] = u_new
            u[end] = t_target
            sol[:, frame] = u
            push!(log, (frame = frame, time = t_target, status = status,
                        iterations = iters, residual = residuals[end], step = dt))
            t = t_target
            popfirst!(queue)
            cutbacks = 0
            continue
        end

        # Не сошлось: дробим подшаг (ANSYS automatic time stepping).
        if cutback >= 1 || cutbacks >= max_cutbacks || dt * cutback < min_step
            @warn "static_solver!: подшаг не пройден" t = t t_target = t_target status = status iterations = iters cutbacks = cutbacks
            break
        end
        cutbacks += 1
        pushfirst!(queue, t + dt * cutback)
    end

    if frame >= 1
        # Хвост — последним кадром, чтобы sol[:, end] всегда был ответом.
        for i in (frame + 1):size(sol, 2)
            sol[:, i] = sol[:, frame]
        end
    else
        for i in 1:size(sol, 2)
            sol[:, i] = u0
        end
        @warn "static_solver!: ни один подшаг не сошёлся, в sol — начальное приближение"
    end

    return log
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