using Pkg; Pkg.activate("./examples")
using Flexia
using ForwardDiff
using GLMakie
using StaticArrays
using CSV
using DataFrames
using GLMakie
using Jurvis
using Statistics
using Peaks
using DSP
using FFTW
using Printf
using LinearAlgebra
# using CSV
# using DataFrames

const OUTDIR     = "Flexia/graphs_22_06_26"


function find_peaks(data::Vector{Float64}; baseline::Float64, threshold::Float64 = 0.5, before::Int = 4000, after::Int = 16000)
    
    window_size = before + after + 1  # всего 301 точка на пик
    peak_indices = Int[]
    
    for i in 1:length(data)
        if data[i] > baseline + threshold &&
           data[i] >= data[i-1] &&
           data[i] >= data[i+1]
            push!(peak_indices, i)
        end
    end
    
    # Убираем слишком близкие пики (берём только самый высокий в окне)
    filtered_peaks = Int[]
    skip_until = 0
    for i in 1:length(peak_indices)
        idx = peak_indices[i]
        if idx <= skip_until
            continue
        end
        # Среди пиков в радиусе `after` берём самый высокий
        group = filter(j -> abs(peak_indices[j] - idx) <= after, 1:length(peak_indices))
        best = peak_indices[group[argmax(data[peak_indices[group]])]]
        push!(filtered_peaks, best)
        skip_until = best + after
    end
    
    # Формируем матрицу: строки = точки окна, столбцы = пики
    n_peaks = length(filtered_peaks)
    result = fill(NaN, window_size, n_peaks)
    
    for i in eachindex(filtered_peaks)
        peak_idx = filtered_peaks[i]
        i_start = peak_idx - before
        i_end   = peak_idx + after
        
        # Определяем, какие строки реально доступны (края массива)
        src_start = max(i_start, 1)
        src_end   = min(i_end, length(data))
        
        row_start = src_start - i_start + 1
        row_end   = row_start + (src_end - src_start)
        
        result[row_start:row_end, i] = data[src_start:src_end]
    end
    
    return result, filtered_peaks
end

function specrum_signal(result::Matrix{Float64})
    fig = Figure(size = (1800, 920))
    xlims!(ax1, 1, 500)
    ylims!(ax1, 0.0001,0.1) 
    rows = length(result[1, :])÷4 +1
    column = 2
    ax = Matrix{Axis}(undef, rows, column)
    for i in 1:rows-1
        for j in 1:column
            ax[i, j] = Axis(fig[i, j], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", xscale = log10, yscale = log10)
        end  
    end
    ax[rows, 1] = Axis(fig[rows, 1:column], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", xscale = log10, yscale = log10)
    ax[rows, 2] = Axis(fig[rows, 1:column], xscale = log10, xticklabelsvisible = false, yticklabelsvisible = false)
    AX = vec(permutedims(ax))
    for k in 1:(size(result, 2)÷2)
        f, T = easyspectrum(result[:, k]; dt = 1.0e-4)
        allpks= findmaxima(T)
        pks, heights = peakheights!(allpks; min=0.007)
        pks, proms = peakproms!(allpks; min = 0.002)
        label_piks = ["$(round(f[pks][i], sigdigits=3)) Hz" for i in 1:length(pks)]
        scatter!(AX[k], f[pks], T[pks], color = :black, markersize = 10, label = "Top freqs")
        label_piks = ["$(round(f[pks][i], sigdigits=3)) Hz" for i in 1:length(pks)]
        text!(AX[k], f[pks], T[pks], text = label_piks, align = (:left, :bottom), offset = (5, 5))
        lines!(AX[k], f, T, color = colors[k])
        xlims!(AX[k], 1, 500)
        ylims!(AX[k], 0.0001,0.1)
        lines!(AX[size(result, 2)÷2+1], f, T, color = colors[k])  
    end
    display(fig)
    return fig, ax
end

function parsing(data::Matrix{Float64},df::DataFrame)
    data[:, 1] = parse.(Float64, string.(replace.(df[!, 1], "," => ".")))
    data[:, 2] = parse.(Float64, string.(replace.(df[!, 2], "," => ".")))
    data[:,  1] = (data[:,  1].- mean(data[:,  1]))./ sensitivity # получаем ускарение
    data[:, 2] = (data[:,  2].- mean(data[:,  2]))./ sensitivity # получаем ускарение
    return data
end

function phase_between(sig_in, sig_out; dt::Real, do_unwrap = false)
    x = collect(float.(sig_in));  x[isnan.(x)] .= 0.0; x .-= mean(x)
    y = collect(float.(sig_out)); y[isnan.(y)] .= 0.0; y .-= mean(y)
    N = min(length(x), length(y))
    X = rfft(x[1:N]); Y = rfft(y[1:N])
    f = rfftfreq(N, 1/dt)
    S12 = conj.(X) .* Y          # взаимный спектр
    φ   = angle.(S12)            # фаза «выход относительно входа»
    amp = abs.(X) .* abs.(Y)     # для маски
    do_unwrap && (φ = DSP.unwrap(φ))
    return f, φ, amp
end

function easyphase(signal::AbstractVector{<:Real}; dt::Real, do_unwrap::Bool = false)
    x = collect(float.(signal))
    x[isnan.(x)] .= 0.0        
    x .-= mean(x)             
    N = length(x)
    X = rfft(x)               
    f = rfftfreq(N, 1/dt)     
    amp = abs.(X)             # модуль ≈ то, что даёт easyspectrum (АЧХ)
    φ   = angle.(X)           # фаза в радианах, диапазон [-π, π]
    do_unwrap && (φ = DSP.unwrap(φ))
    return f, φ, amp
end

function build_6bar_system(hinge_x_offset, K₁, K_hsp)
    g = 9.81

    CURRENT_TIME = 0.
    # K₁ = 0.5
    K₂ = K₁; K₃ = K₂; K₄ = K₂; K₅ = K₂; K₆ = K₁


    bd2_bd6_offset_x = 0.312
    actuators_offset_y = 0.1084
    angle_5 = deg2rad(-110)

    # Массы в кг
    m1 = 24e-3; m2 = 29e-3; m4 = 51e-3; m5 = 52e-3;
    # Моменты инерции в кг·м²
    # 104726.798394 г·мм² = 104726.798394e-9 кг·м²
    I2 = 104726.798394e-9
    # 28187.202983 г·мм² = 28187.202983e-9 кг·м²
    I4 = 36119.84e-9
    F_container = 1
    omega = 0.5
    # Длины в метрах: было 1.5 дм = 0.15 м, 0.75 дм = 0.075 м, 0.5 дм = 0.05 м
    bd1 = Body2D(m5, 0.1e-30, length = 0.15)
    bd2 = Body2D(m2+m1, I2, length = 0.15)
    bd3 = Body2D(m2+m1, I2, length = 0.15)
    bd4 = Body2D(m4+m1, I4, length = 0.1)
    bd5 = Body2D(m2+m1, I2, length = 0.15)
    bd6 = Body2D(m2+m1, I2, length = 0.15)
    bd7 = Body2D(m5, 0.1e-30, length = 0.15)

    slider_rail_1 = Body2D(m2+m1, I2, length = 0.068)
    slider_rail_2 = Body2D(m2+m1, I2, length = 0.068)
    slider_1 = Body2D(m4, I4, length = 0.0255)
    slider_2 = Body2D(m4, I4, length = 0.0255)

    jnt1 = FixedJoint(bd1)
    jnt8 = FixedJoint(bd7)

    jnt2 = HingeJoint(bd1, bd2)
    jnt3 = HingeJoint(bd2, bd3)

    jnt4 = HingeJoint(bd3, bd4)
    jnt5 = HingeJoint(bd4, bd5)

    jnt6 = HingeJoint(bd5, bd6)
    jnt7 = HingeJoint(bd6, bd7)

    # Left side with actuator


    ground_rail_hinge1 = HingeJoint(bd1, slider_rail_1)
    set_position_on_first_body!(ground_rail_hinge1, SA[hinge_x_offset, actuators_offset_y])
    set_position_on_second_body!(ground_rail_hinge1, SA[0., 0.])


    direcrion_SL = SA[-1., 0.]

    x_offset = bd1.length/2 - slider_1.length/2

    slider_ground_slider_1 = SliderJoint(slider_rail_1, slider_1)
    set_position_on_first_body!(slider_ground_slider_1, SA[x_offset, 0.])
    set_position_on_second_body!(slider_ground_slider_1, SA[0., 0.])
    set_direction_on_first_body!(slider_ground_slider_1, direcrion_SL)
    set_direction_on_second_body!(slider_ground_slider_1, direcrion_SL)

    hor_bd2_hinge = HingeJoint(slider_1, bd2)
    set_position_on_first_body!(hor_bd2_hinge, SA[slider_1.length/2, 0])
    set_position_on_second_body!(hor_bd2_hinge, SA[-(bd2.length/2-actuators_offset_y), 0])


    # right side with actuator

    ground_rail_hinge2 = HingeJoint(bd7, slider_rail_2)
    set_position_on_first_body!(ground_rail_hinge2, SA[0., actuators_offset_y])
    set_position_on_second_body!(ground_rail_hinge2, SA[0., 0.])

    hor_bd6_hinge = HingeJoint(bd6, slider_2)
    set_position_on_first_body!(hor_bd6_hinge, SA[-(actuators_offset_y-bd6.length/2), 0])
    set_position_on_second_body!(hor_bd6_hinge, SA[-(slider_2.length/2), 0])

    direcrion_SL_rev = SA[-1., 0.]
    slider_ground_slider_2 = SliderJoint(slider_rail_2, slider_2)
    set_position_on_first_body!(slider_ground_slider_2, SA[-x_offset, 0.])
    set_position_on_second_body!(slider_ground_slider_2, SA[0., 0.])
    set_direction_on_first_body!(slider_ground_slider_2, direcrion_SL_rev)
    set_direction_on_second_body!(slider_ground_slider_2, direcrion_SL_rev)

    tcp1 = TorsionalSpring(jnt2, K₁, deg2rad(-90), 0.1, 0.03)
    tcp2 = TorsionalSpring(jnt3, K₂, deg2rad(45), 0.1, 0.03)

    tcp3 = TorsionalSpring(jnt4, K₃, deg2rad(45), 0.1, 0.03)
    tcp4 = TorsionalSpring(jnt5, K₄, deg2rad(45), 0.1, 0.03)

    tcp5 = TorsionalSpring(jnt6, K₅, deg2rad(45), 0.1, 0.03)
    tcp6 = TorsionalSpring(jnt7, K₆, deg2rad(-90), 0.1, 0.03)

    tcp_ground_rail_1 = TorsionalSpring(ground_rail_hinge1, K₁, deg2rad(0), 0.1, 0.03)
    tcp_ground_rail_2 = TorsionalSpring(ground_rail_hinge2, K₁, deg2rad(0), 0.1, 0.03)

    tcp_hor_bd2 = TorsionalSpring(hor_bd2_hinge, K₁, deg2rad(-90), 0.1, 0.03)
    tcp_hor_bd6 = TorsionalSpring(hor_bd6_hinge, K₁, deg2rad(-90), 0.1, 0.03)

    hsp1 = LinearSpring(slider_ground_slider_1;
                        stiffness = K_hsp,
                        damping = 0.01, vis_r = 0.1, vis_N = 6)
    hsp2 = LinearSpring(slider_ground_slider_2;
                        stiffness = K_hsp, 
                        damping = 0.01, vis_r = 0.1, vis_N = 6)

    # Позиции присоединений: всё было в дм, теперь в м (делим на 10)
    set_position_on_first_body!(jnt2, SA[bd1.length/2, 0.])
    set_position_on_second_body!(jnt2, SA[-bd2.length/2, 0.])

    set_position_on_first_body!(jnt3, SA[bd2.length / 2, 0.])
    set_position_on_second_body!(jnt3, SA[-bd3.length /2, 0.])

    set_position_on_first_body!(jnt4, SA[bd3.length/2, 0])
    set_position_on_second_body!(jnt4, SA[-bd4.length/2, 0])

    set_position_on_first_body!(jnt5, SA[bd4.length/2, 0])
    set_position_on_second_body!(jnt5, SA[-bd5.length/2, 0])

    set_position_on_first_body!(jnt6, SA[bd5.length/2, 0])
    set_position_on_second_body!(jnt6, SA[-bd6.length/2, 0])

    set_position_on_first_body!(jnt7, SA[bd6.length/2, 0])
    set_position_on_second_body!(jnt7, SA[-bd7.length/2, 0])


    # Позиция фиксации тела 7: 6.12 дм = 0.612 м
    jnt8.pos = SA[bd1.length/2 + bd2_bd6_offset_x + bd7.length/2, 0.]

    sys = MBSystem2D()

    add!(sys, bd1)
    add!(sys, bd2)
    add!(sys, bd3)
    add!(sys, bd4)
    add!(sys, bd5)
    add!(sys, bd6)
    add!(sys, bd7)

    add!(sys, slider_rail_1)
    add!(sys, slider_1)

    add!(sys, slider_rail_2)
    add!(sys, slider_2)

    add!(sys, jnt1)

    add!(sys, jnt2)
    add!(sys, jnt3)
    add!(sys, jnt4)
    add!(sys, jnt5)
    add!(sys, jnt6)
    add!(sys, jnt7)

    add!(sys, jnt8)

    add!(sys, ground_rail_hinge1)
    add!(sys, ground_rail_hinge2)

    add!(sys, slider_ground_slider_1)
    add!(sys,slider_ground_slider_2)

    add!(sys,hor_bd2_hinge)
    add!(sys,hor_bd6_hinge)

    add!(sys, tcp1)
    add!(sys, tcp2)
    add!(sys, tcp3)
    add!(sys, tcp4)
    add!(sys, tcp5)
    add!(sys, tcp6)

    add!(sys,tcp_ground_rail_1)
    add!(sys,tcp_ground_rail_2)

    add!(sys,tcp_hor_bd2)
    add!(sys,tcp_hor_bd6)

    add!(sys,hsp1)
    add!(sys,hsp2)

    if (!assemble!(sys))
        println("Assembling failed!")
    end

    func = sys.rhs

    jacoby = (x) -> ForwardDiff.jacobian(func, x)

    bd1_x_ind, bd1_y_ind, bd1_t_ind = get_body_position_dofs(sys, bd1)

    bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(sys, bd2)

    bd3_x_ind, bd3_y_ind, bd3_t_ind = get_body_position_dofs(sys, bd3)

    bd4_x_ind, bd4_y_ind, bd4_t_ind = get_body_position_dofs(sys, bd4)
    bd4_vx_ind, bd4_vy_ind, bd4_vt_ind = get_body_velocity_dofs(sys, bd4)

    bd5_x_ind, bd5_y_ind, bd5_t_ind = get_body_position_dofs(sys, bd5)

    bd6_x_ind, bd6_y_ind, bd6_t_ind = get_body_position_dofs(sys, bd6)

    bd7_x_ind, bd7_y_ind, bd7_t_ind = get_body_position_dofs(sys, bd7)

    bd4_Vx_ind, bd4_Vy_ind, bd4_Vt_ind = get_body_velocity_dofs(sys, bd4)


    initial = zeros(number_of_dofs(sys))


    # Все начальные координаты переведены из дм в м (делением на 10)
    initial[bd1_x_ind] = 0.

    initial[bd2_x_ind] = bd1.length/2
    initial[bd2_y_ind] = bd2.length/2
    initial[bd2_t_ind] = deg2rad(90)

    initial[bd3_x_ind] = bd1.length/2 + bd3.length*cos(deg2rad(45))/2
    initial[bd3_y_ind] = bd2.length + bd3.length*sin(deg2rad(45))/2
    initial[bd3_t_ind] = deg2rad(45)

    initial[bd4_x_ind] = bd1.length/2 + bd3.length*cos(deg2rad(45)) + bd4.length/2
    initial[bd4_y_ind] = bd2.length + bd3.length*sin(deg2rad(45))
    initial[bd4_t_ind] = deg2rad(0)

    initial[bd5_x_ind] = bd1.length/2 + bd2_bd6_offset_x - bd5.length*cos(deg2rad(45))/2
    initial[bd5_y_ind] = bd6.length + bd5.length*sin(deg2rad(45))/2
    initial[bd5_t_ind] = deg2rad(-45)

    initial[bd6_x_ind] = bd1.length/2 + bd2_bd6_offset_x
    initial[bd6_y_ind] = bd6.length/2
    initial[bd6_t_ind] = deg2rad(-90)

    slider_rail_1_X, slider_rail_1_Y, _ =  get_body_position_dofs(sys, slider_rail_1)
    initial[slider_rail_1_X] = 0.
    initial[slider_rail_1_Y] = actuators_offset_y

    slider_rail_2_X, slider_rail_2_Y, _ =  get_body_position_dofs(sys, slider_rail_2)
    initial[slider_rail_2_X] = bd1.length/2 + bd2_bd6_offset_x + bd7.length/2
    initial[slider_rail_2_Y] = actuators_offset_y



    slider_1_X, slider_1_Y, _ =  get_body_position_dofs(sys, slider_1)
    initial[slider_1_X] = bd1.length/2 - slider_1.length/2
    initial[slider_1_Y] = actuators_offset_y

    slider_2_X, slider_2_Y, _ =  get_body_position_dofs(sys, slider_2)
    initial[slider_2_X] = bd1.length/2 + bd2_bd6_offset_x + slider_2.length/2
    initial[slider_2_Y] = actuators_offset_y

    # Начальное время (последний элемент)
    initial[end] = 0.0

    initial[bd7_x_ind] = bd1.length/2 + bd2_bd6_offset_x + bd7.length/2


    f = (t)-> func([initial[1:end-1]; t])[bd4_vx_ind]



    jacoby(initial)
    #

    mass = get_mass_matrix(sys)
    t_end = 10*π / omega
    time_span = LinRange(0,t_end, 501)

    n_steps = length(time_span)

    sol1 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
    cros!(sol1, initial, mass, func, jacoby, step(time_span))
    _, _, bd2_t_ind = get_body_position_dofs(sys, bd2)
    angle = rad2deg.(sol1[bd2_t_ind, end])

    return sys, sol1, initial, angle
end

function find_offset(K₁, K_hsp, target_deg; hinge_x_offset_start = 0, max_iter = 50, tol = 1e-2,  h = 1e-4)
    hinge_x_offset = float(hinge_x_offset_start)

    target_deg_new = 180 - target_deg
    for i in 1:max_iter
        f = build_6bar_system(hinge_x_offset, K₁, K_hsp)[4] - target_deg_new
        if norm(f) < tol
            return hinge_x_offset, "сошлось за $i шагов"
        end
        fp = (build_6bar_system(hinge_x_offset + h, K₁, K_hsp)[4] - build_6bar_system(hinge_x_offset - h, K₁, K_hsp)[4]) / (2h)
        hinge_x_offset = hinge_x_offset - f / fp

    end
    
    return hinge_x_offset, "не сошлось за $max_iter"
end
function apply_stiffnesses!(sys, k, gr)

    for j in eachindex(gr)
        for force in sys.forces[gr[j]]
            force.stiffness = k[j]
        end
    end
    if (!assemble!(sys))
        error("Assembling failed!")
    end

    # вызвать assemble!(sys)
end

function compute_sensitivity_matrix(sys, position, groups)
    initial_freqs = Flexia.find_natural_freqs(sys, position)
    sensitivity_matrix = Matrix{Float64}(undef, length(initial_freqs), length(groups))

    for (i, group) in enumerate(groups)
        for force in sys.forces[group]
            force.stiffness += 1e-6
        end
        assemble!(sys)
        freqs_i = Flexia.find_natural_freqs(sys, position)
        sensitivity_matrix[:, i] = (freqs_i - initial_freqs) / 1e-6
        for force in sys.forces[group]
            force.stiffness -= 1e-6
        end
    end
    assemble!(sys)
    return sensitivity_matrix
end

function update_stiffnesses!( ω_exp, initial_stiffness, stiffness_groups, target_deg; α=0.3, k_min=1e-8, tol=1e-6, max_iter=50)
    hinge_x_offset = 0.
    target_deg_new = 180 - target_deg
    sys, _, position, _ = build_6bar_system(hinge_x_offset, initial_stiffness[1], initial_stiffness[2])
    L = length(stiffness_groups)
    N = length(ω_exp)
    apply_stiffnesses!(sys, initial_stiffness, stiffness_groups)
    new_stiffness = copy(initial_stiffness)
    for m in 1:max_iter

        offset, _ = find_offset(new_stiffness[1], new_stiffness[2], target_deg)
        # println("offset = $offset")
        sys_new, sol_new, _, angle = build_6bar_system(offset, new_stiffness[1], new_stiffness[2])
        q_eq = sol_new[:, end]
        ω_m = Flexia.find_natural_freqs(sys_new, q_eq)[1:N]
        Δω = ω_m - ω_exp
        if norm(Δω) < tol
            @printf "target=%d°: K₁=%.4f, K_hsp=%.2f, ω_m=%s, ω_exp=%s, |Δω|=%.2e\n" target_deg new_stiffness[1] new_stiffness[2] ω_m ω_exp norm(Δω)
            omega = 0.5
            t_end_2 = 10 * π / omega
            time_span = LinRange(0, t_end_2, 501)
            # animate(sys_new, sol_new, time_span, "Flexia/out/new_6bar_test$target_deg.mp4"; framerate = Int(30), limits = (-0.1, 0.6, -0.1, 0.5))
            return new_stiffness, offset, "сошлось за $m шагов"
        end
        T_matrix = compute_sensitivity_matrix(sys_new, q_eq, stiffness_groups)
        T_inv = pinv(T_matrix[:, :])
        new_stiffness -= T_inv[:, 1:N] * Δω

        #ограничеваем шаг
        # for j in 1:L
        #    if abs(Δk[j]) > α * abs(k[j])
        #     Δk[j] = sign(Δk[j]) * α * abs(k[j])
        #    end
        # end
        #что бы не было отрицательных
        for j in 1:L
            new_stiffness[j] = max(new_stiffness[j], k_min)
        end
        apply_stiffnesses!(sys_new, new_stiffness, stiffness_groups)
    end
    return new_stiffness, offset, "не сошлось за $max_iter"
end

function calculate_sistem(F_max, T, t_end, damp_torsional, damp_linear, hinge_x_offset; t_impact = 0.1)
    time_span = LinRange(0.005, t_end, 15001)

    g = 9.81
    CURRENT_TIME = 0.
    K_hsp = 8319.5
    K₁ = 0.16830395
    K₂ = K₁
    K₃ = K₂
    K₄ = K₂
    K₅ = K₂
    K₆ = K₁

    bd2_bd6_offset_x = 0.312
    actuators_offset_y = 0.1084

    # Массы в кг
    m1 = 24e-3
    m2 = 29e-3
    m4 = 51e-3
    m5 = 100000000.

    # Моменты инерции в кг·м²
    I2 = 104726.798394e-9
    I4 = 36119.84e-9
    F_container = 1
    omega = 0.5

    # Длины в метрах
    bd1 = Body2D(m5, 100000000., length = 0.15)
    bd2 = Body2D(m2+m1, I2, length = 0.15)
    bd3 = Body2D(m2+m1, I2, length = 0.15)
    bd4 = Body2D(m4+m1, I4, length = 0.1)
    bd5 = Body2D(m2+m1, I2, length = 0.15)
    bd6 = Body2D(m2+m1, I2, length = 0.15)
    bd7 = Body2D(m5, 100000000., length = 0.15)

    slider_rail_1 = Body2D(m2+m1, I2, length = 0.068)
    slider_rail_2 = Body2D(m2+m1, I2, length = 0.068)

    slider_1 = Body2D(m4, I4, length = 0.0255)
    slider_2 = Body2D(m4, I4, length = 0.0255)

    jnt1 = FixedJoint(bd1)
    jnt8 = FixedJoint(bd7)

    jnt2 = HingeJoint(bd1, bd2)
    jnt3 = HingeJoint(bd2, bd3)

    jnt4 = HingeJoint(bd3, bd4)
    jnt5 = HingeJoint(bd4, bd5)

    jnt6 = HingeJoint(bd5, bd6)
    jnt7 = HingeJoint(bd6, bd7)

    # Left side with actuator
    ground_rail_hinge1 = HingeJoint(bd1, slider_rail_1)
    set_position_on_first_body!(ground_rail_hinge1, SA[hinge_x_offset, actuators_offset_y])
    set_position_on_second_body!(ground_rail_hinge1, SA[0., 0.])

    direcrion_SL = SA[-1., 0.]

    x_offset = bd1.length/2 - slider_1.length/2

    slider_ground_slider_1 = SliderJoint(slider_rail_1, slider_1)
    set_position_on_first_body!(slider_ground_slider_1, SA[x_offset, 0.])
    set_position_on_second_body!(slider_ground_slider_1, SA[0., 0.])
    set_direction_on_first_body!(slider_ground_slider_1, direcrion_SL)
    set_direction_on_second_body!(slider_ground_slider_1, direcrion_SL)

    hor_bd2_hinge = HingeJoint(slider_1, bd2)
    set_position_on_first_body!(hor_bd2_hinge, SA[slider_1.length/2, 0])
    set_position_on_second_body!(hor_bd2_hinge, SA[-(bd2.length/2-actuators_offset_y), 0])

    # Right side with actuator
    ground_rail_hinge2 = HingeJoint(bd7, slider_rail_2)
    set_position_on_first_body!(ground_rail_hinge2, SA[0., actuators_offset_y])
    set_position_on_second_body!(ground_rail_hinge2, SA[0., 0.])

    hor_bd6_hinge = HingeJoint(bd6, slider_2)
    set_position_on_first_body!(hor_bd6_hinge, SA[-(actuators_offset_y-bd6.length/2), 0])
    set_position_on_second_body!(hor_bd6_hinge, SA[-(slider_2.length/2), 0])

    direcrion_SL_rev = SA[-1., 0.]
    slider_ground_slider_2 = SliderJoint(slider_rail_2, slider_2)
    set_position_on_first_body!(slider_ground_slider_2, SA[-x_offset, 0.])
    set_position_on_second_body!(slider_ground_slider_2, SA[0., 0.])
    set_direction_on_first_body!(slider_ground_slider_2, direcrion_SL_rev)
    set_direction_on_second_body!(slider_ground_slider_2, direcrion_SL_rev)

    tcp1 = TorsionalSpring(jnt2, K₁, deg2rad(-90), damp_torsional, 0.03)
    tcp2 = TorsionalSpring(jnt3, K₂, deg2rad(45), damp_torsional, 0.03)

    tcp3 = TorsionalSpring(jnt4, K₃, deg2rad(45), damp_torsional, 0.03)
    tcp4 = TorsionalSpring(jnt5, K₄, deg2rad(45), damp_torsional, 0.03)

    tcp5 = TorsionalSpring(jnt6, K₅, deg2rad(45), damp_torsional, 0.03)
    tcp6 = TorsionalSpring(jnt7, K₆, deg2rad(-90), damp_torsional, 0.03)

    tcp_ground_rail_1 = TorsionalSpring(ground_rail_hinge1, K₁, deg2rad(0), damp_torsional, 0.03)
    tcp_ground_rail_2 = TorsionalSpring(ground_rail_hinge2, K₁, deg2rad(0), damp_torsional, 0.03)

    tcp_hor_bd2 = TorsionalSpring(hor_bd2_hinge, K₁, deg2rad(-90), damp_torsional, 0.03)
    tcp_hor_bd6 = TorsionalSpring(hor_bd6_hinge, K₁, deg2rad(-90), damp_torsional, 0.03)

    hsp1 = LinearSpring(slider_ground_slider_1;
                        stiffness = K_hsp,
                        damping = damp_linear, vis_r = 0.1, vis_N = 6)
    hsp2 = LinearSpring(slider_ground_slider_2;
                        stiffness = K_hsp,
                        damping = damp_linear, vis_r = 0.1, vis_N = 6)

    # Позиции присоединений шарниров
    set_position_on_first_body!(jnt2, SA[bd1.length/2, 0.])
    set_position_on_second_body!(jnt2, SA[-bd2.length/2, 0.])

    set_position_on_first_body!(jnt3, SA[bd2.length / 2, 0.])
    set_position_on_second_body!(jnt3, SA[-bd3.length /2, 0.])

    set_position_on_first_body!(jnt4, SA[bd3.length/2, 0])
    set_position_on_second_body!(jnt4, SA[-bd4.length/2, 0])

    set_position_on_first_body!(jnt5, SA[bd4.length/2, 0])
    set_position_on_second_body!(jnt5, SA[-bd5.length/2, 0])

    set_position_on_first_body!(jnt6, SA[bd5.length/2, 0])
    set_position_on_second_body!(jnt6, SA[-bd6.length/2, 0])

    set_position_on_first_body!(jnt7, SA[bd6.length/2, 0])
    set_position_on_second_body!(jnt7, SA[-bd7.length/2, 0])

    # ── Удар: выключаемый флагом impact_on ──
    # Ref{Bool} позволяет переключать удар после сборки системы и остаётся AD-безопасным.
    impact_on = Ref(false)
    impact_force = (t) -> ifelse(impact_on[] & (t_impact <= t <= t_impact + T),
                                 F_max * sin(pi * (t - t_impact) / T)^2, zero(t))

    hummer_force = Flexia.BodyTimeVariableForce(bd5;
        force = impact_force,
        direction = SA_F64[0, -1],
        position = SA_F64[-0.02, 0.0]
    )

    # Позиция фиксации тела 7
    jnt8.pos = SA[bd1.length/2 + bd2_bd6_offset_x + bd7.length/2, 0.]

    eola_x  = FrameLocalAccelerationSensor(body=bd4, position=SA_F64[0.015, 0.0],  axis=:x)
    eola2_y = FrameLocalAccelerationSensor(body=bd4, position=SA_F64[-0.015, 0.0], axis=:y)

    sys = MBSystem2D()

    add!(sys, bd1)
    add!(sys, bd2)
    add!(sys, bd3)
    add!(sys, bd4)
    add!(sys, bd5)
    add!(sys, bd6)
    add!(sys, bd7)

    add!(sys, slider_rail_1)
    add!(sys, slider_1)

    add!(sys, slider_rail_2)
    add!(sys, slider_2)

    add!(sys, jnt1)

    add!(sys, jnt2)
    add!(sys, jnt3)
    add!(sys, jnt4)
    add!(sys, jnt5)
    add!(sys, jnt6)
    add!(sys, jnt7)

    add!(sys, jnt8)

    add!(sys, ground_rail_hinge1)
    add!(sys, ground_rail_hinge2)

    add!(sys, slider_ground_slider_1)
    add!(sys, slider_ground_slider_2)

    add!(sys, hor_bd2_hinge)
    add!(sys, hor_bd6_hinge)

    add!(sys, tcp1)
    add!(sys, tcp2)
    add!(sys, tcp3)
    add!(sys, tcp4)
    add!(sys, tcp5)
    add!(sys, tcp6)

    add!(sys, tcp_ground_rail_1)
    add!(sys, tcp_ground_rail_2)

    add!(sys, tcp_hor_bd2)
    add!(sys, tcp_hor_bd6)

    add!(sys, hsp1)
    add!(sys, hsp2)
    add!(sys, hummer_force)

    add!(sys, eola_x)
    add!(sys, eola2_y)

    if (!assemble!(sys))
        println("Assembling failed!")
    end

    func   = sys.rhs
    jacoby = (x) -> ForwardDiff.jacobian(func, x)

    bd1_x_ind, bd1_y_ind, bd1_t_ind = get_body_position_dofs(sys, bd1)
    bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(sys, bd2)
    bd3_x_ind, bd3_y_ind, bd3_t_ind = get_body_position_dofs(sys, bd3)
    bd4_x_ind, bd4_y_ind, bd4_t_ind = get_body_position_dofs(sys, bd4)
    bd4_vx_ind, bd4_vy_ind, bd4_vt_ind = get_body_velocity_dofs(sys, bd4)
    bd5_x_ind, bd5_y_ind, bd5_t_ind = get_body_position_dofs(sys, bd5)
    bd6_x_ind, bd6_y_ind, bd6_t_ind = get_body_position_dofs(sys, bd6)
    bd7_x_ind, bd7_y_ind, bd7_t_ind = get_body_position_dofs(sys, bd7)

    initial = zeros(number_of_dofs(sys))

    # Начальные координаты (прямая геометрия — стартовая догадка для релаксации)
    initial[bd1_x_ind] = 0.

    initial[bd2_x_ind] = bd1.length/2
    initial[bd2_y_ind] = bd2.length/2
    initial[bd2_t_ind] = deg2rad(90)

    initial[bd3_x_ind] = bd1.length/2 + bd3.length*cos(deg2rad(45))/2
    initial[bd3_y_ind] = bd2.length + bd3.length*sin(deg2rad(45))/2
    initial[bd3_t_ind] = deg2rad(45)

    initial[bd4_x_ind] = bd1.length/2 + bd3.length*cos(deg2rad(45)) + bd4.length/2
    initial[bd4_y_ind] = bd2.length + bd3.length*sin(deg2rad(45))
    initial[bd4_t_ind] = deg2rad(0)

    initial[bd5_x_ind] = bd1.length/2 + bd2_bd6_offset_x - bd5.length*cos(deg2rad(45))/2
    initial[bd5_y_ind] = bd6.length + bd5.length*sin(deg2rad(45))/2
    initial[bd5_t_ind] = deg2rad(-45)

    initial[bd6_x_ind] = bd1.length/2 + bd2_bd6_offset_x
    initial[bd6_y_ind] = bd6.length/2
    initial[bd6_t_ind] = deg2rad(-90)

    slider_rail_1_X, slider_rail_1_Y, _ = get_body_position_dofs(sys, slider_rail_1)
    initial[slider_rail_1_X] = 0.
    initial[slider_rail_1_Y] = actuators_offset_y

    slider_rail_2_X, slider_rail_2_Y, _ = get_body_position_dofs(sys, slider_rail_2)
    initial[slider_rail_2_X] = bd1.length/2 + bd2_bd6_offset_x + bd7.length/2
    initial[slider_rail_2_Y] = actuators_offset_y

    slider_1_X, slider_1_Y, _ = get_body_position_dofs(sys, slider_1)
    initial[slider_1_X] = bd1.length/2 - slider_1.length/2
    initial[slider_1_Y] = actuators_offset_y

    slider_2_X, slider_2_Y, _ = get_body_position_dofs(sys, slider_2)
    initial[slider_2_X] = bd1.length/2 + bd2_bd6_offset_x + slider_2.length/2
    initial[slider_2_Y] = actuators_offset_y

    initial[end] = 0.0   # начальное время
    initial[bd7_x_ind] = bd1.length/2 + bd2_bd6_offset_x + bd7.length/2

    mass = get_mass_matrix(sys)

    # ── Фаза 1: релаксация в равновесие (удар выключен) ──
    # Долгий грубый прогон CROS численно гасит быстрые моды и осаживает
    # механизм в подогнутую позу (заданную hinge_x_offset).
    impact_on[] = false
    relax_span = LinRange(0.0, 10π / omega, 501)
    sol_relax = Matrix{Float64}(undef, number_of_dofs(sys), length(relax_span))
    cros!(sol_relax, initial, mass, func, jacoby, step(relax_span))

    q_eq = copy(sol_relax[:, end])   # равновесие: подогнутая поза, скорости ≈ 0
    q_eq[end] = 0.0                  # сброс времени к началу фазы удара

    # ── Фаза 2: удар из равновесия ──
    # Старт из покоя в равновесии → до t_impact baseline плоский,
    # вся динамика порождается только молотком.
    impact_on[] = true
    sol = simulate(sys, q_eq, time_span)

    sol1 = sol[:, 3:end]
    time_span_sp = LinRange(0, t_end, length(time_span)-2)
    meas = similar(sol1, (2, length(time_span_sp)))
    buf = meas[:, 1]
    for i in axes(meas, 2)
        sys.measure!(buf, sol1[:, i])
        meas[:, i] .= buf
    end

    _, _, bd2_t_ind = get_body_position_dofs(sys, bd2)
    angle = rad2deg.(sol1[bd2_t_ind, end])

    return sys, sol1, time_span, time_span_sp, meas, angle
end

begin # парсинг данных и поиск пиков
    
    colors = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :yellow, :black, :gray]
    VD = 5.0 # напр питания
    sensitivity = 0.020 # чувствительность В/g (20 мВ/g)
    df_80=CSV.read("Flexia/examples/23.05.26/new-80-ok.csv", DataFrame)
    df_90=CSV.read("Flexia/examples/23.05.26/new-90-ok.csv", DataFrame)
    df_100=CSV.read("Flexia/examples/23.05.26/new-100-ok.csv", DataFrame)
    t_80 = 1:length(df_80[:, 1])
    t_90 = 1:length(df_90[:, 1])
    t_100 = 1:length(df_100[:, 1])
    dt = 1.0e-4
    time_80 = t_80.* dt
    time_90 = t_90.* dt
    time_100 = t_100.* dt

    signal_80 = Matrix{Float64}(undef, length(df_80[:, 1]), 2)
    signal_90 = Matrix{Float64}(undef, length(df_90[:, 1]), 2)
    signal_100 = Matrix{Float64}(undef, length(df_100[:, 1]), 2)

    data_80 = parsing(signal_80, df_80)
    data_90 = parsing(signal_90, df_90)
    data_100 = parsing(signal_100, df_100)

    result_80, peaks_80 = find_peaks(data_80[:,1], baseline = data_80[1,1])
    println("Найдено пиков: ", length(peaks_80))
    result_80 = hcat(result_80, find_peaks(data_80[:,2], baseline = data_80[1,2])[1][:, :]) # result_80 - матрица где первые 4 столбца это пики для первого акселерометра, а с 5 по 8 для второго акселерометра
    peaks_80 =hcat(peaks_80, find_peaks(data_80[:,2], baseline = data_80[1,2])[2][:, :])

    # удар 2 и 4 не полные
    result_90, peaks_90 = find_peaks(data_90[:,1], baseline = data_90[1,1])
    println("Найдено пиков: ", length(peaks_90))
    result_90 = hcat(result_90, find_peaks(data_90[:,2], baseline = data_90[1,2])[1][:, :]) # result_90 - матрица где первые 2 столбца это пики для первого акселерометра, а с 3 по 4 для второго акселерометра
    peaks_90 =hcat(peaks_90, find_peaks(data_90[:,2], baseline = data_90[1,2])[2][:, :])


    # удар 2 и 4 не полные
    result_100, peaks_100 = find_peaks(data_100[:,1], baseline = data_100[1,1])
    println("Найдено пиков: ", length(peaks_100))
    result_100 = hcat(result_100, find_peaks(data_100[:,2], baseline = data_100[1,2])[1][:, :]) # result_100 - матрица где первые 2 столбца это пики для первого акселерометра, а с 3 по 4 для второго акселерометра
    peaks_100 =hcat(peaks_100, find_peaks(data_100[:,2], baseline = data_100[1,2])[2][:, :])
end
##

damp_torsional = 1.0e-3
damp_linear = 1.0e-4
F_max = 20.8
T = 0.007
t_end = 2.1

# Расчет жесткости и сдвига пружины для каждой конфигурации (80°, 90°, 100°) на основе экспериментальных частот.
# k = [0.01, 5000.0]
# ω_exp_80 = [8, 33]
# ω_exp_90 = [7, 30]# эксп первая и вторая частоты
# ω_exp_100 = [6, 36]
# gr1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
# gr2 = [11, 12]
# gr = [gr1, gr2]

# k_result_80, offset_80, status_80 = update_stiffnesses!( ω_exp_80, k, gr, 80; α=0.3, max_iter=4500, tol=1e-3)
# k_result_90, offset_90, status_90 = update_stiffnesses!( ω_exp_90, k, gr, 90; α=0.3, max_iter=4500, tol=1e-3)
# k_result_100, offset_100, status_100 = update_stiffnesses!( ω_exp_100, k, gr, 100; α=0.3, max_iter=4500, tol=1e-3)

hinge_x_offset_80 = -0.0189329 # для 80°
hinge_x_offset_90 = 0.0 # для 90°
hinge_x_offset_100 = 0.01902 # для 100°

sys, sol1, time_span, time_span_sp, meas, angle = calculate_sistem(F_max, T, t_end, damp_torsional, damp_linear, hinge_x_offset_80)
signal_x_model = meas[1, :]
signal_y_model = meas[2, :]
signal_x_model[1:9] .= 0.0;
signal_y_model[1:8] .= 0.0;

animate(sys, sol1, time_span_sp, "./new_6bar.mp4"; framerate = 1 ÷ (1*step(time_span)), limits = (-0.1, 0.6, -0.1, 0.5))

# animate(sys, sol1[:, 1:150], time_span_sp[1:150], "./new_6bar_slow.mp4"; framerate = 1 ÷ (40*step(time_span)), limits = (-0.1, 0.6, -0.1, 0.5))
# Flexia.draw_static(sys, initial; limits = (-0.1, 0.6, -0.1, 0.5))

# fig = Figure()
# ax = Axis(fig[1, 1])
# l1 = lines!(ax, time_span_sp, -meas[1,:], linestyle = :dashdot)
# l2 = lines!(ax, time_span_sp, meas[2,:])
# Legend(fig[1 , 2], [l1,l2], ["A1X", "A1Y"], framevisible = false, halign = :right, valign = :top)
# fig


signal_x_experiment = result_90[:, 7]
signal_y_experiment = result_90[:, 4]

Flexia.find_natural_freqs(sys, sol1[:, end])
begin #графики спектра симуляции и удара
    ω = [7,30, 35.78]
    dt_model = step(time_span_sp)
    fx_model, Tx_model = easyspectrum(signal_x_model; dt = dt_model)
    fy_model, Ty_model = easyspectrum(signal_y_model; dt = dt_model)

    fx_exp, Tx_exp = easyspectrum(signal_x_experiment; dt = dt)
    fy_exp, Ty_exp = easyspectrum(signal_y_experiment; dt = dt)

    fig = Figure(size = (800, 600))
    ax1 = Axis(fig[1,1], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", title = "Signal_X", yscale = log10)
    ax2 = Axis(fig[2,1], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", title = "Signal_Y", yscale = log10)
    l1 = lines!(ax1, fx_exp, Tx_exp, color=:red)
    lines!(ax2, fy_exp, Ty_exp, color=:red)

    l2 = lines!(ax1, fx_model, Tx_model, color=:blue)
    lines!(ax2, fy_model, Ty_model, color=:blue)

    l3 = vlines!(ax1, ω, color=:black, linestyle = :dashdot)
    vlines!(ax2, ω, color=:black, linestyle = :dashdot)
    xlims!(ax1, 0, 100)
    xlims!(ax2, 0, 100)
    ylims!(ax1, 0.0001, 1)
    ylims!(ax2, 0.0001, 1)
    Legend(fig[1 , 2], [l1, l2, l3], ["experiment", "model", "experimental modes"], framevisible = false, halign = :right, valign = :top)
    fig
end
##

fig = Figure(size = (800, 600))
ax1 = Axis(fig[1,1], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_Y")
lines!(ax1, ((1:length(signal_x_experiment)).*dt).- 0.38, signal_x_experiment, color=:green, linestyle = :dash)
xlims!(ax1, 0, 0.5)
ylims!(ax1, -5, 5)
fig
begin #графики симуляции и удара

    fig = Figure(size = (800, 600))
    ax1 = Axis(fig[1,1], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_Y")
    ax2 = Axis(fig[2,1], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_X")

    lines!(ax1, ((1:length(signal_x_experiment)).*dt).- 0.38, signal_x_experiment, color=:green, linestyle = :dash)
    lines!(ax2, ((1:length(signal_y_experiment)).*dt).- 0.3975, signal_y_experiment, color=:green, linestyle = :dash)

    lines!(ax1, time_span_sp .-0.005, signal_x_model, color=:red)
    lines!(ax2, time_span_sp .-0.005, signal_y_model, color=:red)

    xlims!(ax1, 0, 0.5)
    ylims!(ax1, -5, 5)
    xlims!(ax2, 0, 0.5)
    ylims!(ax2, -5, 5)

    fig
end


# Flexia.find_natural_freqs(sys, initial)



begin
    fig = Figure()
    ax1 = Axis(fig[1,1], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_11_Y")
    ax2 = Axis(fig[2,1], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_12_Y")
    ax3 = Axis(fig[3,1], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_13_Y")
    ax4 = Axis(fig[4,1], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_14_Y")

    lines!(ax1, (1:length(result_90[:, 1])).*dt, result_90[:, 1], color=:red)
    lines!(ax2, (1:length(result_90[:, 2])).*dt, result_90[:, 2], color=:blue)
    lines!(ax3, (1:length(result_90[:, 3])).*dt, result_90[:, 3], color=:green)
    lines!(ax4, (1:length(result_90[:, 4])).*dt, result_90[:, 4], color=:orange)

    ax5 = Axis(fig[1,2], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_11_X")
    ax6 = Axis(fig[2,2], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_12_X")
    ax7 = Axis(fig[3,2], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_13_X")
    ax8 = Axis(fig[4,2], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_14_X")

    lines!(ax5, (1:length(result_90[:, 5])).*dt, result_90[:, 5], color=:red)
    lines!(ax6, (1:length(result_90[:, 6])).*dt, result_90[:, 6], color=:blue)
    lines!(ax7, (1:length(result_90[:, 7])).*dt, result_90[:, 7], color=:green)
    lines!(ax8, (1:length(result_90[:, 8])).*dt, result_90[:, 8], color=:orange)
    fig
end




fig = Figure()
ax1 = Axis(fig[1,1], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_1_Y")
ax2 = Axis(fig[2,1], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_1_Y")
lines!(ax1, (1:length(result_90[:, 1])).*dt, result_90[:, 1], color=:red)
lines!(ax2, (1:length(result_90[:, 2])).*dt, result_90[:, 5], color=:blue)



