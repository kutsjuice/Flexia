using Pkg; Pkg.activate("./examples")
using Flexia
using ForwardDiff
using GLMakie
using StaticArrays
using CSV
using DataFrames
const g = 9.81

CURRENT_TIME = 0.

# K₁ = 0.03
# K₂ = 0.0115
K₁ = 0.02
K₂ = 0.02
K₃ = K₂
K₄ = K₂
K₅ = K₂
K₆ = K₁


const angle_5 = deg2rad(-110)

# Массы в кг
m1 = 23e-3
m2 = 29e-3
m4 = 15e-3

m5 = 52e-3
# Моменты инерции в кг·м²
# 104726.798394 г·мм² = 104726.798394e-9 кг·м²
I2 = 104726.798394e-9
# 28187.202983 г·мм² = 28187.202983e-9 кг·м²
I4 = 28187.202983e-9
F_container = 20.
omega = 2. 
# Длины в метрах: было 1.5 дм = 0.15 м, 0.75 дм = 0.075 м, 0.5 дм = 0.05 м
bd1 = Body2D(m5, 0.1e-30, length = 0.15)
bd2 = Body2D(m2+m1, I2, length = 0.075)
bd3 = Body2D(m2+m1, I2, length = 0.075)

bd4 = Body2D(m4+m1, I4, length = 0.05)
# bd4.forces[1] = (x) -> -(m4+m1) * F_container
bd4.forces[1] = (x) -> 0.05 * sin(2.0 * x[7])
bd5 = Body2D(m2+m1, I2, length = 0.075)
bd6 = Body2D(m2+m1, I2, length = 0.075)
bd7 = Body2D(m5, 0.1e-30, length = 0.15)

jnt1 = FixedJoint(bd1)
jnt8 = FixedJoint(bd7)

jnt2 = HingeJoint(bd1, bd2)
jnt3 = HingeJoint(bd2, bd3)

jnt4 = HingeJoint(bd3, bd4)
jnt5 = HingeJoint(bd4, bd5)

jnt6 = HingeJoint(bd5, bd6)
jnt7 = HingeJoint(bd6, bd7)

tcp1 = TorsionalSpring(bd1, bd2, K₁, deg2rad(-90), 0.)
tcp2 = TorsionalSpring(bd2, bd3, K₂, deg2rad(45), 0.)

tcp3 = TorsionalSpring(bd3, bd4, K₃, deg2rad(45), 0.)
tcp4 = TorsionalSpring(bd4, bd5, K₄, deg2rad(45), 0.)

tcp5 = TorsionalSpring(bd5, bd6, K₅, deg2rad(45), 0.)
tcp6 = TorsionalSpring(bd6, bd7, K₆, angle_5, 0.)

# Позиции присоединений: всё было в дм, теперь в м (делим на 10)
set_position_on_first_body!(jnt2, SA[0.15, 0.])
set_position_on_second_body!(jnt2, SA[-0.075, 0.])

set_position_on_first_body!(jnt3, SA[0.075, 0])
set_position_on_second_body!(jnt3, SA[-0.075, 0])

set_position_on_first_body!(jnt4, SA[0.075, 0])
set_position_on_second_body!(jnt4, SA[-0.05, 0])

set_position_on_first_body!(jnt5, SA[0.05, 0])
set_position_on_second_body!(jnt5, SA[-0.075, 0])

set_position_on_first_body!(jnt6, SA[0.075, 0])
set_position_on_second_body!(jnt6, SA[-0.075, 0])

set_position_on_first_body!(jnt7, SA[0.075, 0])
set_position_on_second_body!(jnt7, SA[-0.15, 0])


# Позиция фиксации тела 7: 6.12 дм = 0.612 м
jnt8.pos = SA[0.612, 0.]

sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)
add!(sys, bd3)
add!(sys, bd4)
add!(sys, bd5)
add!(sys, bd6)
add!(sys, bd7)

add!(sys, jnt1)

add!(sys, jnt2)
add!(sys, jnt3)
add!(sys, jnt4)
add!(sys, jnt5)
add!(sys, jnt6)
add!(sys, jnt7)

add!(sys, jnt8)

add!(sys, tcp1)
add!(sys, tcp2)
add!(sys, tcp3)
add!(sys, tcp4)
add!(sys, tcp5)
add!(sys, tcp6)

add_time!(sys)

if (!assemble!(sys))
    println("Assembling failed!")
end

func = sys.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)

bd1_x_ind, bd1_y_ind, bd1_t_ind = get_body_position_dofs(sys, bd1)

bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(sys, bd2)

bd3_x_ind, bd3_y_ind, bd3_t_ind = get_body_position_dofs(sys, bd3)

bd4_x_ind, bd4_y_ind, bd4_t_ind = get_body_position_dofs(sys, bd4)

bd5_x_ind, bd5_y_ind, bd5_t_ind = get_body_position_dofs(sys, bd5)

bd6_x_ind, bd6_y_ind, bd6_t_ind = get_body_position_dofs(sys, bd6)

bd7_x_ind, bd7_y_ind, bd7_t_ind = get_body_position_dofs(sys, bd7)

bd4_Vx_ind, bd4_Vy_ind, bd4_Vt_ind = get_body_velocity_dofs(sys, bd4)


initial = zeros(number_of_dofs(sys))


# Все начальные координаты переведены из дм в м (делением на 10)
initial[bd1_x_ind] = 0.

initial[bd2_x_ind] = 0.15
initial[bd2_y_ind] = 0.15/2
initial[bd2_t_ind] = deg2rad(90)

initial[bd3_x_ind] = 0.15 + 0.15*cos(deg2rad(45))/2
initial[bd3_y_ind] = 0.15 + 0.15*sin(deg2rad(45))/2
initial[bd3_t_ind] = deg2rad(45)

initial[bd4_x_ind] = 0.15 + 0.05 + 0.15*cos(deg2rad(45))
initial[bd4_y_ind] = 0.15 + 0.15*sin(deg2rad(45))
initial[bd4_t_ind] = deg2rad(0)

initial[bd5_x_ind] = 0.462 - 0.15*cos(deg2rad(46))/2
initial[bd5_y_ind] = 0.15 + 0.15*sin(deg2rad(45))/2
initial[bd5_t_ind] = deg2rad(-45)

initial[bd6_x_ind] = 0.462
initial[bd6_y_ind] = 0.15/2
initial[bd6_t_ind] = deg2rad(-90)

# Начальное время (последний элемент)
initial[end] = 0.0

# # === ЯВНЫЙ ИНТЕГРАТОР RK4 ===
# dt = 0.0005
# time_span = 0:dt:10.0
# n_steps = length(time_span)

# sol = Matrix{Float64}(undef, number_of_dofs(sys), n_steps)
# sol[:, 1] = initial

# # Классический RK4
# for i in 2:n_steps
#     u = sol[:, i-1]
    
#     k1 = func(u)
#     k2 = func(u + 0.5*dt*k1)
#     k3 = func(u + 0.5*dt*k2)
#     k4 = func(u + dt*k3)
    
#     sol[:, i] = u + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
# end

# animate(sys, sol, time_span, "bd4_harmonic_force.mp4"; 
#     framerate = 30, limits = (-5, 5, -5, 5))

# println("Симуляция завершена! Файл: bd4_harmonic_force.mp4")
initial[bd7_x_ind] = 0.612

func(initial)
jacoby(initial)

mass = zeros(number_of_dofs(sys), number_of_dofs(sys));
for i in 1:last_body_dof(sys)
    mass[i, i] = 1
end
mass[end, end] = 1.0
time_span = range(0., 20., 200)
length(time_span)
sol_euler = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
sol_euler[:,1] = initial
n_steps = length(time_span)
length(time_span)
sol1 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
cros!(sol1, initial, mass, func, jacoby, step(time_span))

# Лимиты графика тоже в метрах
animate(sys, sol1, time_span, "triv_dyn.mp4"; framerate = 30, limits = (0., 0.7, 0., 0.7))

# # sol2 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
# for i in 2:n_steps
#     global  CURRENT_TIME = time_span[i-1]
#     rhs_val = func(sol_euler[:, i-1])
#     sol_euler[:, i] = sol_euler[:, i-1] + dt * rhs_val
# end

# # static_solver!(sol2, initial, func, jacoby)
display(sol1[bd4_Vx_ind,:])
# display(sol1[bd4_y_ind,:])
# # animate(sys, sol2, time_span, "triv_stat.mp4"; framerate = 60, limits = (-0.1, 0.5, -0.1, 0.5))

f1 = Figure()
ax1 = Axis(f1[1,1], xscale = identity, aspect = AxisAspect(2),title="Траектория X",xlabel="Сила, Н",ylabel="Перемещение, мм")
ylims!(minimum(sol1[bd4_x_ind,:]), maximum(sol1[bd4_x_ind,:]))
xlims!(minimum(sol1[bd4_Vx_ind,:]), maximum(sol1[bd4_Vx_ind,:]))
l1 = lines!(ax1, sol1[bd4_Vx_ind,:], sol1[bd4_x_ind,:], linestyle = :dot)
save("рисунок траектории динамика флексия.png",f1)
# # # 
recurdyn_data = CSV.File("volna.csv") |> DataFrame

# f2 = Figure()
# ax2 = Axis(f2[1,1], xscale = identity, aspect = AxisAspect(2),title="Траектория X",xlabel="Сила, Н",ylabel="Перемещение, мм")
# ylims!(minimum(recurdyn_data[1:48,2]), maximum(recurdyn_data[1:48,2]))
# xlims!(minimum(recurdyn_data[1:48,3]), maximum(recurdyn_data[1:48,3]))
# l2 = lines!(ax2, recurdyn_data[1:48,3], recurdyn_data[1:48,2], linestyle = :dot)
# save("рисунок траектории динамика рекурдин.png",f2)
perem = [0.]
rec_perem = [0.]
for i in 2:200
    push!(perem, (sol1[bd4_x_ind, i] - sol1[bd4_x_ind,1])*571.4378712)
    push!(rec_perem, (recurdyn_data[i,2] - recurdyn_data[1,2]))
end
perem
rec_perem
 maximum(rec_perem)
  maximum(perem)
f3 = Figure()
ax3 = Axis(f3[1,1], xscale = identity, aspect = AxisAspect(2),title="Перемещение от времени",xlabel="Время, С",ylabel="Перемещение, мм")
ylims!(minimum(perem), maximum(perem))
xlims!(minimum(time_span), maximum(time_span))
l3 = lines!(ax3, time_span, perem)
l4 = lines!(ax3, time_span, rec_perem)
save("перемещение от времени.png", f3)