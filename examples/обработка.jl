using CSV
using DataFrames
using GLMakie
using Jurvis
using Statistics
using Peaks
using DSP
using FFTW
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
const OUTDIR     = "graphs_23_05_06"

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

# графики сигналов
# fig = Figure()
# ax = Axis(fig[1,1], xlabel = "time, [s]", ylabel = "acceleration, [g]")
# lines!(ax, time_80, data_80[:, 2], color=:red, label="acc1")
# fig = Figure()
# ax = Axis(fig[1,1], xlabel = "time, [s]", ylabel = "acceleration, [g]")
# lines!(ax, time_90, data_90[:, 2], color=:blue, label="acc1")
# fig = Figure()
# ax = Axis(fig[1,1], xlabel = "time, [s]", ylabel = "acceleration, [g]")
# lines!(ax, time_100, data_100[:, 1], color=:green, label="acc1")


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

# график удара 


#фильтрация шума
# noise_level = data_80[1:60000, 1]
# time_noise = time[1:60000]

# f_noise, A_noise = easyspectrum(noise_level; dt = 1.0e-4)
# fig = Figure()
# fs = 1/dt  # частота дискретизации (1/dt)
# f1, T1 = easyspectrum(result_80[:, 1]; dt = 1.0e-4)
# lowpass = digitalfilter(Lowpass(50.0), Butterworth(4); fs=fs)
# impact_signal = result_80[:, 1]  # или другой столбец удара
# impact_filtered = filtfilt(lowpass, impact_signal)
# filtered_data = filtfilt(lowpass, data_80[:, 1])
# filtered_frequency, filtered_amplitude = easyspectrum(impact_filtered; dt = 1.0e-4)

# fig = Figure()
# ax1 = Axis(fig[1,1], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", title = "Спектр шума", xscale = log10, yscale = log10)
# ax2 = Axis(fig[1,2], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", title = "без фильтра", xscale = log10, yscale = log10)
# ax3 = Axis(fig[1,3], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", title = "с фильтром", xscale = log10, yscale = log10)
# lines!(ax1, f_noise, A_noise, color=:red)
# lines!(ax2, f1, T1, color=:red)
# lines!(ax2, filtered_frequency, filtered_amplitude, color=:blue)
# lines!(ax3, filtered_frequency, filtered_amplitude, color=:blue)
# xlims!(ax2, 1, 500)
# ylims!(ax2, 0.0001,0.02)


# specrum_signal(result_80)
let
    f_2, φ_2, amp_2 = easyphase(result_80[:, 6]; dt = 1.0e-4, do_unwrap = true)
    f_1, φ_1, amp_1 = easyphase(result_80[:, 2]; dt = 1.0e-4, do_unwrap = true)
    # фаза имеет смысл только там, где амплитуда заметна;
    # на «пустых» частотах это просто шум — маскируем
    thr_1 = 0.02 * maximum(amp_1)
    thr_2 = 0.02 * maximum(amp_2)
    φ_masked_1 = [a > thr_1 ? p : NaN for (p, a) in zip(φ_1, amp_1)]
    φ_masked_2 = [a > thr_2 ? p : NaN for (p, a) in zip(φ_2, amp_2)]
    rng_1 = f_1 .> 0              # частота 0 не дружит с log-осью
    rng_2 = f_2 .> 0              # частота 0 не дружит с log-осью

    fig = Figure(size = (1000, 600))
    ax1 = Axis(fig[1, 1],
        xlabel = "Frequency, [Hz]", ylabel = "Phase, [rad]",
        title = "ФЧХ, 80°, signal 1", xscale = log10)
    ax2 = Axis(fig[1, 2],
        xlabel = "Frequency, [Hz]", ylabel = "Phase, [rad]",
        title = "ФЧХ, 80°, signal 2", xscale = log10)
    lines!(ax2, f_2[rng_2], φ_2[rng_2], color = :blue)
    xlims!(ax2, 1, 500)
    ax1.yticks = (-π:π/2:π, ["-π", "-π/2", "0", "π/2", "π"])
    lines!(ax1, f_1[rng_1], φ_1[rng_1], color = :red)
    xlims!(ax1, 1, 500)
    ax1.yticks = (-π:π/2:π, ["-π", "-π/2", "0", "π/2", "π"])
    display(fig)
    ylims!(ax1, -4, 4)
    ylims!(ax2, -4, 4)
end
let
    onset = 4000   # before из find_peaks: пик стоит на отсчёте onset+1

    # обрезаем окно к моменту удара (убираем линейный набег) + считаем фазу с маской
    function phase_of(col)
        sig = result_80[onset+1:end, col]
        sig = sig[.!isnan.(sig)]                       # на случай неполного удара
        f, φ, amp = easyphase(sig; dt = 1.0e-4, do_unwrap = true)
        thr = 0.02 * maximum(amp)
        φm = [(a > thr && fi > 0) ? p : NaN            # fi > 0: частота 0 не дружит с log
              for (p, a, fi) in zip(φ, amp, f)]
        return f, φm
    end

    f1, φ1 = phase_of(2)   # signal 1, удар 2
    f2, φ2 = phase_of(6)   # signal 2, удар 2

    fig = Figure(size = (1000, 600))
    ax1 = Axis(fig[1, 1], xlabel = "Frequency, [Hz]", ylabel = "Phase, [rad]",
               title = "ФЧХ, 80°, signal 1", xscale = log10)
    ax2 = Axis(fig[1, 2], xlabel = "Frequency, [Hz]", ylabel = "Phase, [rad]",
               title = "ФЧХ, 80°, signal 2", xscale = log10)

    lines!(ax1, f1, φ1, color = :red)
    lines!(ax2, f2, φ2, color = :blue)

    for ax in (ax1, ax2)
        xlims!(ax, 1, 500)
        ylims!(ax, -10, 10)
    end
    linkyaxes!(ax1, ax2)        # одинаковый масштаб Y → графики сравнимы

    display(fig)
end


begin # спектр для удара для 80 градусов датчик 1 и 2.
    fig = Figure(size = (1800, 920))
    ax1 = Axis(fig[1,1], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", title = "ch_1_soft2_80°")
    ax2 = Axis(fig[1,2], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", title = "ch_2_soft2_80°")
    f2, T2 = easyspectrum(result_80[:, 2]; dt = 1.0e-4)
    f6, T6 = easyspectrum(result_80[:, 6]; dt = 1.0e-4)
    lines!(ax1, f2, T2/T2[2], color = colors[2])
    lines!(ax2, f6, T6/T6[2], color = colors[3])
    xlims!(ax1, 5, 35)
    ylims!(ax1, 1,35)  
    xlims!(ax2, 5, 35)
    ylims!(ax2, 1, 50)
    display(fig)   
end
 

let# спектр для удара для 80 градусов датчик 1. выбрал 3 удар

    fig = Figure(size = (1800, 920))
    ax1 = Axis(fig[1,1], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", xscale = log10, yscale = log10)
    ax2 = Axis(fig[1,2], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", xscale = log10, yscale = log10,)
    ax3 = Axis(fig[2,1], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", xscale = log10, yscale = log10)
    ax4 = Axis(fig[2,2], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", xscale = log10, yscale = log10)
    ax5 = Axis(fig[3,1:2], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", xscale = log10, yscale = log10)

    f1, T1 = easyspectrum(result_80[:, 1]; dt = 1.0e-4)
    lines!(ax1, f1, T1, color = colors[1])
    xlims!(ax1, 1, 500)
    ylims!(ax1, 0.0001,0.1)   
    allpks1= findmaxima(T1)
    pks1, heights1 = peakheights!(allpks1; min=0.007)
    pks1, proms1 = peakproms!(allpks1; min = 0.002)
    scatter!(ax1, f1[pks1], T1[pks1], color = :black, markersize = 10, label = "Top freqs")
    label_piks = ["$(round(f1[pks1][i], sigdigits=3)) Hz" for i in 1:length(pks1)]
    text!(ax1, f1[pks1], T1[pks1], text = label_piks, align = (:left, :bottom), offset = (5, 5))


    f2, T2 = easyspectrum(result_80[:, 2]; dt = 1.0e-4)
    lines!(ax2, f2, T2, color = colors[2])
    xlims!(ax2, 1, 500)
    ylims!(ax2, 0.0001,0.1)  
    allpks2= findmaxima(T2)
    pks2, heights = peakheights!(allpks2; min=0.01)
    pks2, proms = peakproms!(allpks2; min = 0.002)

    scatter!(ax2, f2[pks2], T2[pks2], color = :black, markersize = 10, label = "Top freqs")
    label_piks = ["$(round(f2[pks2][i], sigdigits=3)) Hz" for i in 1:length(pks2)]
    text!(ax2, f2[pks2], T2[pks2], text = label_piks, align = (:left, :bottom), offset = (5, 5))


    # plotpeaks!(plt, f2, T2, removed; show_prominences=false, show_widths=false, markercolor=:Gray) # removed peaks

    f3, T3 = easyspectrum(result_80[:, 3]; dt = 1.0e-4)
    lines!(ax3, f3, T3, color = colors[3])
    xlims!(ax3, 1, 500)
    ylims!(ax3, 0.0001,0.1)  
    allpks3= findmaxima(T3)
    pks3, heights = peakheights!(allpks3; min=0.0092)
    pks3, proms = peakproms!(allpks3; min = 0.002)
    # peaksplot(f3, T3, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax3, f3[pks3], T3[pks3], color = :black, markersize = 10, label = "Top freqs")
    label_piks = ["$(round(f3[pks3][i], sigdigits=3)) Hz" for i in 1:length(pks3)]
    text!(ax3, f3[pks3], T3[pks3], text = label_piks, align = (:left, :bottom), offset = (5, 5))

    f4, T4 = easyspectrum(result_80[:, 4]; dt = 1.0e-4)
    lines!(ax4, f4, T4, color = colors[4])
    xlims!(ax4, 1,500)   # подберите диапазон под ваш сигнал
    ylims!(ax4, 0.0001,0.1)
    allpks4= findmaxima(T4)
    pks4, heights = peakheights!(allpks4; min=0.0112)
    pks4, proms = peakproms!(allpks4; min = 0.002)
    # peaksplot(f4, T4, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax4, f4[pks4], T4[pks4], color = :black, markersize = 10, label = "Top freqs")
    label_piks = ["$(round(f4[pks4][i], sigdigits=3)) Hz" for i in 1:length(pks4)]
    text!(ax4, f4[pks4], T4[pks4], text = label_piks, align = (:left, :bottom), offset = (5, 5))

    lines!(ax5, f1, T1, color = colors[1])
    lines!(ax5, f2, T2, color = colors[2])
    lines!(ax5, f3, T3, color = colors[3])
    lines!(ax5, f4, T4, color = colors[4])

    fname = joinpath(OUTDIR, "summary_all_spectr_80_1.png")
    save(fname, fig)

    display(fig)
end

let # спектр для удара для 80 градусов датчик 2. выбрал 3 удар
    fig = Figure(size = (1800, 920))
    ax1 = Axis(fig[1,1], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", xscale = log10, yscale = log10)
    ax2 = Axis(fig[1,2], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", xscale = log10, yscale = log10,)
    ax3 = Axis(fig[2,1], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", xscale = log10, yscale = log10)
    ax4 = Axis(fig[2,2], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", xscale = log10, yscale = log10)
    ax5 = Axis(fig[3,1:2], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", xscale = log10, yscale = log10)

    f1, T1 = easyspectrum(result_80[:, 5]; dt = 1.0e-4)
    lines!(ax1, f1, T1, color = colors[1])
    xlims!(ax1, 1, 500)
    ylims!(ax1, 0.0001,0.1)   
    allpks1= findmaxima(T1)
    pks1, heights1 = peakheights!(allpks1; min=0.007)
    pks1, proms1 = peakproms!(allpks1; min = 0.002)
    scatter!(ax1, f1[pks1], T1[pks1], color = :black, markersize = 10, label = "Top freqs")
    label_piks = ["$(round(f1[pks1][i], sigdigits=3)) Hz" for i in 1:length(pks1)]
    text!(ax1, f1[pks1], T1[pks1], text = label_piks, align = (:left, :bottom), offset = (5, 5))


    f2, T2 = easyspectrum(result_80[:, 6]; dt = 1.0e-4)
    lines!(ax2, f2, T2, color = colors[2])
    xlims!(ax2, 1, 500)
    ylims!(ax2, 0.0001,0.1)  
    allpks2= findmaxima(T2)
    pks2, heights = peakheights!(allpks2; min=0.01)
    pks2, proms = peakproms!(allpks2; min = 0.002)

    scatter!(ax2, f2[pks2], T2[pks2], color = :black, markersize = 10, label = "Top freqs")
    label_piks = ["$(round(f2[pks2][i], sigdigits=3)) Hz" for i in 1:length(pks2)]
    text!(ax2, f2[pks2], T2[pks2], text = label_piks, align = (:left, :bottom), offset = (5, 5))


    # plotpeaks!(plt, f2, T2, removed; show_prominences=false, show_widths=false, markercolor=:Gray) # removed peaks

    f3, T3 = easyspectrum(result_80[:, 7]; dt = 1.0e-4)
    lines!(ax3, f3, T3, color = colors[3])
    xlims!(ax3, 1, 500)
    ylims!(ax3, 0.0001,0.1)  
    allpks3= findmaxima(T3)
    pks3, heights = peakheights!(allpks3; min=0.0092)
    pks3, proms = peakproms!(allpks3; min = 0.002)
    # peaksplot(f3, T3, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax3, f3[pks3], T3[pks3], color = :black, markersize = 10, label = "Top freqs")
    label_piks = ["$(round(f3[pks3][i], sigdigits=3)) Hz" for i in 1:length(pks3)]
    text!(ax3, f3[pks3], T3[pks3], text = label_piks, align = (:left, :bottom), offset = (5, 5))

    f4, T4 = easyspectrum(result_80[:, 8]; dt = 1.0e-4)
    lines!(ax4, f4, T4, color = colors[4])
    xlims!(ax4, 1,500)   # подберите диапазон под ваш сигнал
    ylims!(ax4, 0.0001,0.1)
    allpks4= findmaxima(T4)
    pks4, heights = peakheights!(allpks4; min=0.0112)
    pks4, proms = peakproms!(allpks4; min = 0.002)
    # peaksplot(f4, T4, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax4, f4[pks4], T4[pks4], color = :black, markersize = 10, label = "Top freqs")
    label_piks = ["$(round(f4[pks4][i], sigdigits=3)) Hz" for i in 1:length(pks4)]
    text!(ax4, f4[pks4], T4[pks4], text = label_piks, align = (:left, :bottom), offset = (5, 5))

    lines!(ax5, f1, T1, color = colors[1])
    lines!(ax5, f2, T2, color = colors[2])
    lines!(ax5, f3, T3, color = colors[3])
    lines!(ax5, f4, T4, color = colors[4])

    fname = joinpath(OUTDIR, "summary_all_spectr_80_2.png")
    save(fname, fig)

    display(fig)
end

let # спектр для удара 80
    fig = Figure(size = (1800, 420))
    ax1 = Axis(fig[1,1], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", xscale = log10, yscale = log10, title = "Signal 1")
    ax2 = Axis(fig[1,2], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", xscale = log10, yscale = log10, title = "Signal 2")
    Label(fig[0, :], "angle 80°", fontsize = 22, font = :bold)
    f1, T1 = easyspectrum(result_80[:, 3]; dt = 1.0e-4)
    lines!(ax1, f1, T1, color = colors[1])
    xlims!(ax1, 1, 500)
    ylims!(ax1, 0.0001,0.1)  
    allpks1= findmaxima(T1)
    pks1, heights = peakheights!(allpks1; min=0.0092)
    pks1, proms = peakproms!(allpks1; min = 0.002)
    # peaksplot(f1, T1, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax1, f1[pks1], T1[pks1], color = :black, markersize = 10, label = "Top freqs")
    label_piks = ["$(round(f1[pks1][i], sigdigits=3)) Hz" for i in 1:length(pks1)]
    text!(ax1, f1[pks1], T1[pks1], text = label_piks, align = (:left, :bottom), offset = (5, 5))
   
    f2, T2 = easyspectrum(result_80[:, 7]; dt = 1.0e-4)
    lines!(ax2, f2, T2, color = colors[2])
    xlims!(ax2, 1, 500)
    ylims!(ax2, 0.0001,0.1)  
    allpks2= findmaxima(T2)
    pks2, heights = peakheights!(allpks2; min=0.0112)
    pks2, proms = peakproms!(allpks2; min = 0.002)
    # peaksplot(f2, T2, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax2, f2[pks2], T2[pks2], color = :black, markersize = 10, label = "Top freqs")
    label_piks = ["$(round(f2[pks2][i], sigdigits=3)) Hz" for i in 1:length(pks2)]
    text!(ax2, f2[pks2], T2[pks2], text = label_piks, align = (:left, :bottom), offset = (5, 5))
    fname = joinpath(OUTDIR, "spectr_80.png")
    save(fname, fig)
    display(fig)

end

let # спектр для удара для 80 градусов оба датчика вместе
    fig = Figure(size = (900, 900))
    ax1 = Axis(fig[1,1], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", title = "80°", xscale = log10, yscale = log10)
    # ax2= Axis(fig[1,2], xlabel = "Frequency, [Hz]", ylabel = "Amplitude",title = "Signal 2", xscale = log10, yscale = log10)
    # Label(fig[0, 1], "angle 80°", fontsize = 22, font = :bold)
    f1, T1 = easyspectrum(result_80[:, 3]; dt = 1.0e-4)
    lines!(ax1, f1, T1, color = colors[1], label = "signal 1")
    xlims!(ax1, 1,600)   # подберите диапазон под ваш сигнал
    ylims!(ax1, 0.0001,0.1)
    allpks1= findmaxima(T1)
    pks1, heights = peakheights!(allpks1; min=0.0112)
    pks1, proms = peakproms!(allpks1; min = 0.002)
    # peaksplot(f1, T1, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax1, marker = :cross,  f1[pks1], T1[pks1], color = :black, markersize = 10, label = "Top freqs of signal 1")
    label_piks = ["$(round(f1[pks1][i], sigdigits=3)) Hz" for i in 1:length(pks1)]
    text!(ax1, f1[pks1], T1[pks1], text = label_piks, align = (:left, :bottom), offset = (5, 5))
    
    f2, T2 = easyspectrum(result_80[:, 7]; dt = 1.0e-4)
    lines!(ax1,linestyle = :dash, f2, T2, color = colors[2], label = "signal 2")
    allpks2= findmaxima(T2)
    pks2, heights = peakheights!(allpks2; min=0.0005)
    pks2, proms = peakproms!(allpks2; min = 0.0017)
    # peaksplot(f1, T1, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax1, f2[pks2], T2[pks2], color = :black, markersize = 10, label = "Top freqs of signal 2")
    label_piks = ["$(round(f2[pks2][i], sigdigits=3)) Hz" for i in 1:length(pks2)]
    text!(ax1, f2[pks2], T2[pks2], text = label_piks, align = (:left, :bottom), offset = (5, 5))
    axislegend(ax1)
    display(fig)
    fname = joinpath(OUTDIR, "spectr_80_combined.png")
    save(fname, fig)
end

let # спектр для удара для 90 градусов оба датчика отдельно
    fig = Figure(size = (1800, 420))
    ax1 = Axis(fig[1,1], xlabel = "Frequency, [Hz]", ylabel = "Amplitude",title = "Signal 1", xscale = log10, yscale = log10)
    ax2= Axis(fig[1,2], xlabel = "Frequency, [Hz]", ylabel = "Amplitude",title = "Signal 2", xscale = log10, yscale = log10)
    Label(fig[0, :], "angle 90°", fontsize = 22, font = :bold)
    f1, T1 = easyspectrum(result_90[:, 1]; dt = 1.0e-4)
    lines!(ax1, f1, T1, color = colors[4])
    xlims!(ax1, 1,600)   # подберите диапазон под ваш сигнал
    ylims!(ax1, 0.0001,0.1)
    allpks1= findmaxima(T1)
    pks1, heights = peakheights!(allpks1; min=0.0112)
    pks1, proms = peakproms!(allpks1; min = 0.002)
    # peaksplot(f1, T1, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax1, f1[pks1], T1[pks1], color = :black, markersize = 10, label = "Top freqs")
    label_piks = ["$(round(f1[pks1][i], sigdigits=3)) Hz" for i in 1:length(pks1)]
    text!(ax1, f1[pks1], T1[pks1], text = label_piks, align = (:left, :bottom), offset = (5, 5))
    
    f2, T2 = easyspectrum(result_90[:, 3]; dt = 1.0e-4)
    lines!(ax2, f2, T2, color = colors[4])
    xlims!(ax2, 1,600)   # подберите диапазон под ваш сигнал
    ylims!(ax2, 0.0001,0.1)
    allpks2= findmaxima(T2)
    pks2, heights = peakheights!(allpks2; min=0.0005)
    pks2, proms = peakproms!(allpks2; min = 0.0017)
    # peaksplot(f1, T1, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax2, f2[pks2], T2[pks2], color = :black, markersize = 10, label = "Top freqs")
    label_piks = ["$(round(f2[pks2][i], sigdigits=3)) Hz" for i in 1:length(pks2)]
    text!(ax2, f2[pks2], T2[pks2], text = label_piks, align = (:left, :bottom), offset = (5, 5))
    display(fig)
    fname = joinpath(OUTDIR, "spectr_90.png")
    save(fname, fig)
end

let # спектр для удара для 90 градусов оба датчика вместе
    fig = Figure(size = (900, 900))
    ax1 = Axis(fig[1,1], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", title = "90°", xscale = log10, yscale = log10)
    # ax2= Axis(fig[1,2], xlabel = "Frequency, [Hz]", ylabel = "Amplitude",title = "Signal 2", xscale = log10, yscale = log10)
    # Label(fig[0, 1], "angle 90°", fontsize = 22, font = :bold)
    f1, T1 = easyspectrum(result_90[:, 1]; dt = 1.0e-4)
    lines!(ax1, f1, T1, color = colors[1], label = "signal 1")
    xlims!(ax1, 1,600)   # подберите диапазон под ваш сигнал
    ylims!(ax1, 0.0001,0.1)
    allpks1= findmaxima(T1)
    pks1, heights = peakheights!(allpks1; min=0.0112)
    pks1, proms = peakproms!(allpks1; min = 0.002)
    # peaksplot(f1, T1, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax1, marker = :cross,  f1[pks1], T1[pks1], color = :black, markersize = 10, label = "Top freqs of signal 1")
    label_piks = ["$(round(f1[pks1][i], sigdigits=3)) Hz" for i in 1:length(pks1)]
    text!(ax1, f1[pks1], T1[pks1], text = label_piks, align = (:left, :bottom), offset = (5, 5))
    
    f2, T2 = easyspectrum(result_90[:, 3]; dt = 1.0e-4)
    lines!(ax1,linestyle = :dash, f2, T2, color = colors[2], label = "signal 2")
    allpks2= findmaxima(T2)
    pks2, heights = peakheights!(allpks2; min=0.0005)
    pks2, proms = peakproms!(allpks2; min = 0.0017)
    # peaksplot(f1, T1, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax1, f2[pks2], T2[pks2], color = :black, markersize = 10, label = "Top freqs of signal 2")
    label_piks = ["$(round(f2[pks2][i], sigdigits=3)) Hz" for i in 1:length(pks2)]
    text!(ax1, f2[pks2], T2[pks2], text = label_piks, align = (:left, :bottom), offset = (5, 5))
    axislegend(ax1)
    display(fig)
    fname = joinpath(OUTDIR, "spectr_90_combined.png")
    save(fname, fig)
end

let # спектр для удара для 110 градусов оба датчика отдельно
    fig = Figure(size = (1800, 420))
    ax1 = Axis(fig[1,1], xlabel = "Frequency, [Hz]", ylabel = "Amplitude",title = "Signal 1", xscale = log10, yscale = log10)
    ax2= Axis(fig[1,2], xlabel = "Frequency, [Hz]", ylabel = "Amplitude",title = "Signal 2", xscale = log10, yscale = log10)
    Label(fig[0, :], "angle 110°", fontsize = 22, font = :bold)
    f1, T1 = easyspectrum(result_110[:, 1]; dt = 1.0e-4)
    lines!(ax1, f1, T1, color = colors[4])
    xlims!(ax1, 1,500)   # подберите диапазон под ваш сигнал
    ylims!(ax1, 0.0001,0.1)
    allpks1= findmaxima(T1)
    pks1, heights = peakheights!(allpks1; min=0.005)
    pks1, proms = peakproms!(allpks1; min = 0.002)
    # peaksplot(f1, T1, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax1, f1[pks1], T1[pks1], color = :black, markersize = 10, label = "Top freqs")
    label_piks = ["$(round(f1[pks1][i], sigdigits=3)) Hz" for i in 1:length(pks1)]
    text!(ax1, f1[pks1], T1[pks1], text = label_piks, align = (:left, :bottom), offset = (5, 5))
    
    f2, T2 = easyspectrum(result_110[:, 3]; dt = 1.0e-4)
    lines!(ax2, f2, T2, color = colors[4])
    xlims!(ax2, 1,500)   # подберите диапазон под ваш сигнал
    ylims!(ax2, 0.0001,0.1)
    allpks2= findmaxima(T2)
    pks2, heights = peakheights!(allpks2; min=0.005)
    pks2, proms = peakproms!(allpks2; min = 0.002)
    # peaksplot(f1, T1, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax2, f2[pks2], T2[pks2], color = :black, markersize = 10, label = "Top freqs")
    label_piks = ["$(round(f2[pks2][i], sigdigits=3)) Hz" for i in 1:length(pks2)]
    text!(ax2, f2[pks2], T2[pks2], text = label_piks, align = (:left, :bottom), offset = (5, 5))
    display(fig)
    fname = joinpath(OUTDIR, "spectr_110.png")
    save(fname, fig)
end

let # спектр для удара для 110 градусов оба датчика вместе
    fig = Figure(size = (900, 900))
    ax1 = Axis(fig[1,1], xlabel = "Frequency, [Hz]", ylabel = "Amplitude", title = "110°", xscale = log10, yscale = log10)
    # ax2= Axis(fig[1,2], xlabel = "Frequency, [Hz]", ylabel = "Amplitude",title = "Signal 2", xscale = log10, yscale = log10)
    # Label(fig[0, 1], "angle 110°", fontsize = 22, font = :bold)
    f1, T1 = easyspectrum(result_110[:, 1]; dt = 1.0e-4)
    lines!(ax1, f1, T1, color = colors[1], label = "signal 1")
    xlims!(ax1, 1,600)   # подберите диапазон под ваш сигнал
    ylims!(ax1, 0.0001,0.1)
    allpks1= findmaxima(T1)
    pks1, heights = peakheights!(allpks1; min=0.0112)
    pks1, proms = peakproms!(allpks1; min = 0.002)
    # peaksplot(f1, T1, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax1, marker = :cross,  f1[pks1], T1[pks1], color = :black, markersize = 10, label = "Top freqs of signal 1")
    label_piks = ["$(round(f1[pks1][i], sigdigits=3)) Hz" for i in 1:length(pks1)]
    text!(ax1, f1[pks1], T1[pks1], text = label_piks, align = (:left, :bottom), offset = (5, 5))
    
    f2, T2 = easyspectrum(result_110[:, 3]; dt = 1.0e-4)
    lines!(ax1,linestyle = :dash, f2, T2, color = colors[2], label = "signal 2")
    allpks2= findmaxima(T2)
    pks2, heights = peakheights!(allpks2; min=0.0005)
    pks2, proms = peakproms!(allpks2; min = 0.0017)
    # peaksplot(f1, T1, significant_peaks; show_prominences=true, show_widths=false, color=:red, markersize=8)
    # подберите диапазон под ваш сигнал
    scatter!(ax1, f2[pks2], T2[pks2], color = :black, markersize = 10, label = "Top freqs of signal 2")
    label_piks = ["$(round(f2[pks2][i], sigdigits=3)) Hz" for i in 1:length(pks2)]
    text!(ax1, f2[pks2], T2[pks2], text = label_piks, align = (:left, :bottom), offset = (5, 5))
    axislegend(ax1)
    display(fig)
    fname = joinpath(OUTDIR, "spectr_110_combined.png")
    save(fname, fig)
end

let# графики ударов
    fig1 = Figure(size = (1800, 420))
    ax1 = Axis(fig1[1,1], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_1")
    ax2 = Axis(fig1[1,2], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_2")
    t_80 = 1:length(result_80[:, 3])
    t_80 = t_80 .* dt
    Label(fig1[0, :], "angle 80°", fontsize = 22, font = :bold)
    lines!(ax1, t_80, result_80[:, 3], color=:blue)
    xlims!(ax1, 0.39,0.7)
    lines!(ax2, t_80, result_80[:, 7], color=:blue)
    xlims!(ax2, 0.39,0.7)
    fname = joinpath(OUTDIR, "punch_80.png")
    save(fname, fig1)

    t_90 = 1:length(result_90[:, 3])
    t_90 = t_90 .* dt
    fig2 = Figure(size = (1800, 320))
    ax1 = Axis(fig2[1,1], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_1")
    ax2 = Axis(fig2[1,2], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_2")
    Label(fig2[0, :], "angle 90°", fontsize = 22, font = :bold) 
    lines!(ax1, t_90, result_90[:, 1], color=:blue)
    xlims!(ax1, 0.39,0.7)
    lines!(ax2, t_90, result_90[:, 3], color=:blue)
    xlims!(ax2, 0.39,0.7)
    fname = joinpath(OUTDIR, "punch_90.png")
    save(fname, fig2)

    t_110 = 1:length(result_110[:, 3])
    t_110 = t_110 .* dt
    fig3 = Figure(size = (1800, 320))
    ax1 = Axis(fig3[1,1], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_1")
    ax2 = Axis(fig3[1,2], xlabel = "time, [s]", ylabel = "acceleration, [g]", title = "Signal_2")
    Label(fig3[0, :], "angle 110°", fontsize = 22, font = :bold)
    lines!(ax1, t_110, result_110[:, 1], color=:blue)
    xlims!(ax1, 0.39,0.7)
    lines!(ax2, t_110, result_110[:, 3], color=:blue)
    xlims!(ax2, 0.39,0.7)
    fname = joinpath(OUTDIR, "punch_110.png")
    save(fname, fig3)
    GLMakie.closeall()
    display(GLMakie.Screen(), fig1)
    display(GLMakie.Screen(), fig2)
    display(GLMakie.Screen(), fig3)
end

function plotchannel(data::MeasuredData; ch::Int64 = 1)
    fig = Figure()
    ax  = Axis(fig[1, 1],
               xlabel = data.xdata_name,
               ylabel = data.ydata_names[ch],
               title  = data.name)
    lines!(ax, data.xdata, data.ydata[:, ch], label = data.ydata_names[ch])
    axislegend(ax)
    return fig
end

function plotall(data::MeasuredData; skip = Int[], leg = false)
    fig = Figure()
    ax  = Axis(fig[1, 1],
               xlabel = data.xdata_name,
               title  = data.name)
    for i in 1:data.num
        i in skip && continue
        lines!(ax, data.xdata, data.ydata[:, i], label = data.ydata_names[i])
    end
    leg && axislegend(ax)
    display(fig)
    return fig
end

function plotrelation(data::MeasuredData, chX::Int, chY::Int)
    fig = Figure()
    ax  = Axis(fig[1, 1],
               xlabel = data.ydata_names[chX],
               ylabel = data.ydata_names[chY],
               title  = "Dependency of $(data.ydata_names[chY]) from $(data.ydata_names[chX])")
    lines!(ax, data.ydata[:, chX], data.ydata[:, chY])
    return fig
end

time_impact = (1:size(result_80, 1)).* dt
SIGNAL = MeasuredData("impacts_80", time_impact, result_80, "time, [s]", ["acc1_imp1", "acc1_imp2", "acc1_imp3","acc1_imp4", "acc2_imp1", "acc2_imp2", "acc2_imp3", "acc2_imp4"], size(result_80, 1), size(result_80, 2), fill("", size(result_80, 2)))
f = 8
SSA_win = trunc(Int, 5 / (dt * f))
decomposed_data = decomposeSSA(SIGNAL, 1, 10, SSA_win)
spectrum1 = spectrumall(decomposed_data);
spectrum1

# number = 10
# fig = Figure(size = (900, 950))
# ax = Axis(fig[1, 1], xlabel = "time", xscale = log10, yscale = log10)
# lines!(ax, spectrum1.xdata, spectrum1.ydata[:, number], color=colors[number])
# xlims!(ax, 1, 1000)
# ylims!(ax, 0.0001,30)

groups_8hz = [[4, 5]]
grouped_modes_8hz = groupmodes(decomposed_data, groups_8hz)
residual_signal = grouped_modes_8hz.ydata[:, 1] - grouped_modes_8hz.ydata[:, 2]
f_res, T_res = easyspectrum(residual_signal; dt = dt)

let
fig_residual = Figure(size = (1000, 800))
ax_time = Axis(fig_residual[1, 1], 
               xlabel = "time, [s]", 
               ylabel = "acceleration, [g]", 
               title = "Original signal vs signal without 8 Hz mode")
lines!(ax_time, time_impact, grouped_modes_8hz.ydata[:, 1], color = (:gray, 0.5), label = "Original signal")
lines!(ax_time, time_impact, residual_signal, color = :blue, label = "Signal without 8 Hz", linestyle = :dash)
axislegend(ax_time)
ax_spec = Axis(fig_residual[2, 1], 
               xlabel = "Frequency, [Hz]", 
               ylabel = "Amplitude", 
               title = "Spectrum of the signal after subtracting the 8 Hz mode", 
               xscale = log10, 
               yscale = log10)

lines!(ax_spec, f_res, T_res, color = :blue)
lines!(ax_spec, spectrum1.xdata, spectrum1.ydata[:, 1], color = :red)
display(fig_residual)
end
# fig = Figure(size = (900, 950))
# ax = Axis(fig[1, 1], xlabel = "time", )
# lines!(ax, grouped_modes_8hz.xdata, grouped_modes_8hz.ydata[:,3])
# lines!(ax, decomposed_data.xdata, decomposed_data.ydata[:,9])

