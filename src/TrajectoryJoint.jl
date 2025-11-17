abstract type AbstractTrajectoryJoint <: AbstractJoint2D end

mutable struct TrajectoryJoint <: AbstractTrajectoryJoint
    body::Body2D
    trajectory::Function  # функция траектории (t) -> [x, y, θ]
    start_time::Float64   # время начала движения
    duration::Float64     # продолжительность движения
    index::Int64
    
    function TrajectoryJoint(body::Body2D, trajectory::Function, start_time::Float64=0.0, duration::Float64=1.0)
        return new(body, trajectory, start_time, duration, -1)
    end
end

number_of_dofs(::TrajectoryJoint) = 3

function add!(sys::MBSystem2D, joint::TrajectoryJoint)
    push!(sys.joints, joint)
    last_joint_dof = last_lm_dof(sys) + number_of_dofs(joint)
    push!(sys.lmdofs, last_joint_dof)
    setid!(joint, length(sys.joints))
end

function get_lms(sys::MBSystem2D, joint::TrajectoryJoint)
    last_lm = sys.lmdofs[joint.index]
    return SA[
        last_lm-2,
        last_lm-1,
        last_lm-0,
    ]
end

function add_joint_to_rhs!(rhs, state, sys::MBSystem2D, joint::TrajectoryJoint)
    body = joint.body
    last_body_dof = sys.bodiesdofs[body.index]
    
    position_dofs = [
        last_body_dof - 5,
        last_body_dof - 4,
        last_body_dof - 3,
    ]
    
    velocity_dofs = [
        last_body_dof - 2,
        last_body_dof - 1,
        last_body_dof,
    ]
    
    joint_dofs = get_lms(sys, joint)
    
    # Получаем текущее время из состояния (предполагаем, что время хранится в последнем элементе)
    current_time = state[end]  # или использовать глобальное время
    
    if current_time >= joint.start_time && current_time <= joint.start_time + joint.duration
        # Вычисляем желаемую позицию и скорость
        t = current_time - joint.start_time
        desired_pos = joint.trajectory(t)
        
        # Численное дифференцирование для скорости
        dt = 1e-6
        if t + dt <= joint.duration
            pos_plus = joint.trajectory(t + dt)
            desired_vel = (pos_plus - desired_pos) / dt
        else
            pos_minus = joint.trajectory(t - dt)
            desired_vel = (desired_pos - pos_minus) / dt
        end
        
        # Добавляем управляющие моменты
        rhs[velocity_dofs[1]] += state[joint_dofs[1]]
        rhs[velocity_dofs[2]] += state[joint_dofs[2]]
        rhs[velocity_dofs[3]] += state[joint_dofs[3]]
        
        # Ограничения для следования траектории
        rhs[joint_dofs[1]] = state[position_dofs[1]] - desired_pos[1]
        rhs[joint_dofs[2]] = state[position_dofs[2]] - desired_pos[2]
        rhs[joint_dofs[3]] = state[position_dofs[3]] - desired_pos[3]
    else
        # Вне временного интервала - свободное движение
        rhs[velocity_dofs[1]] += state[joint_dofs[1]]
        rhs[velocity_dofs[2]] += state[joint_dofs[2]]
        rhs[velocity_dofs[3]] += state[joint_dofs[3]]
        
        rhs[joint_dofs[1]] = 0.0
        rhs[joint_dofs[2]] = 0.0
        rhs[joint_dofs[3]] = 0.0
    end
end

# Вспомогательные функции для создания траекторий
function circular_trajectory(center, radius, angular_velocity)
    return (t) -> [
        center[1] + radius * cos(angular_velocity * t),
        center[2] + radius * sin(angular_velocity * t),
        angular_velocity * t
    ]
end

function linear_trajectory(start_pos, end_pos, duration)
    return (t) -> [
        start_pos[1] + (end_pos[1] - start_pos[1]) * t / duration,
        start_pos[2] + (end_pos[2] - start_pos[2]) * t / duration,
        start_pos[3] + (end_pos[3] - start_pos[3]) * t / duration
    ]
end

function sinusoidal_trajectory(base_pos, amplitude, frequency, axis=1)
    return (t) -> [
        base_pos[1] + (axis == 1 ? amplitude * sin(frequency * t) : 0.0),
        base_pos[2] + (axis == 2 ? amplitude * sin(frequency * t) : 0.0),
        base_pos[3] + (axis == 3 ? amplitude * sin(frequency * t) : 0.0)
    ]
end

function polynomial_trajectory(coefficients)
    return (t) -> [
        sum(coefficients[1][i] * t^(i-1) for i in 1:length(coefficients[1])),
        sum(coefficients[2][i] * t^(i-1) for i in 1:length(coefficients[2])),
        sum(coefficients[3][i] * t^(i-1) for i in 1:length(coefficients[3]))
    ]
end