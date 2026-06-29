
function measure(state::Vector{Float64}, dstate::Vector{Float64}, sys::MBSystem2D, sensor::AbstractSensor2D)::Float64
    return NaN64
end

function add!(sys::MBSystem2D, sensor::S) where S <: AbstractSensor2D
    push!(sys.sensors, sensor)
    return nothing
end

mutable struct FramePositionSensor{A <: AbstractAxis} <: AbstractPositionSensor
    body::Body2D
    position::SVector{2, Float64}
    axis::A
end
function FramePositionSensor(;body::Body2D, position::SVector{2, Float64}, axis::Union{Int64, Symbol})
    if axis == 1 || :x || axis == :X
        return FramePositionSensor(body, position, XAxis())
    elseif axis == 2 || :y || axis == :Y
        return FramePositionSensor(body, position, YAxis())
    else
        throw(ArgumentError("Axis must be either 1 (for X) or  2 (for Y) "))
    end
end

function measure(state::Vector{Float64}, dstate::Vector{Float64}, sys::MBSystem2D, sensor::FramePositionSensor{XAxis})::Float64
    body = sensor.body;
    dofs = get_body_position_dofs(sys, body);
    x_, y_, r_ = state[dofs];
    return first(SA_F64[cos(r_) -sin(r_)] * sensor.position) + x_;
end


function measure(state::Vector{Float64}, dstate::Vector{Float64}, sys::MBSystem2D, sensor::FramePositionSensor{YAxis})::Float64
    body = sensor.body;
    dofs = get_body_position_dofs(sys, body);
    x_, y_, r_ = state[dofs];
    return first(SA_F64[sin(r_) cos(r_)] * sensor.position) + y_;
end
mutable struct FrameGlobalVelocitySensor{A <: AbstractAxis} <: AbstractSensor2D
    body::Body2D
    position::SVector{2, Float64}
    axis::A
end
function FrameGlobalVelocitySensor(;body::Body2D, position::SVector{2, Float64}, axis::Union{Int64, Symbol})
    if axis == 1 || :x || axis == :X
        return FrameGlobalVelocitySensor(body, position, XAxis())
    elseif axis == 2 || :y || axis == :Y
        return FrameGlobalVelocitySensor(body, position, YAxis())
    else
        throw(ArgumentError("Axis must be either 1 (for X) or  2 (for Y) "))
    end
end

function measure(state::Vector{Float64}, dstate::Vector{Float64}, sys::MBSystem2D, sensor::FrameGlobalVelocitySensor{XAxis})::Float64
    body = sensor.body;
    pos_dofs = get_body_position_dofs(sys, body);
    vel_dofs = get_body_velocity_dofs(sys, body);
    x_, y_, r_ = state[pos_dofs];
    vx_, vy_, vr_ = state[vel_dofs];
    p_x, p_y = sensor.position;

    return vx_ - vr_ * (sin(r_) * p_x + cos(r_) * p_y);
end
function measure(state::Vector{Float64}, dstate::Vector{Float64}, sys::MBSystem2D, sensor::FrameGlobalVelocitySensor{YAxis})::Float64
    body = sensor.body;
    pos_dofs = get_body_position_dofs(sys, body);
    vel_dofs = get_body_velocity_dofs(sys, body);
    x_, y_, r_ = state[pos_dofs];
    vx_, vy_, vr_ = state[vel_dofs];
    p_x, p_y = sensor.position;
    return vy_ + vr_ * (cos(r_) * p_x - sin(r_) * p_y);
end

function add!(sys::MBSystem2D, sensor::FrameGlobalVelocitySensor)
    push!(sys.sensors, sensor)
    return nothing
end


mutable struct FrameLocalVelocitySensor{A<:AbstractAxis} <: AbstractSensor2D
    body::Body2D
    position::SVector{2, Float64}
    rot::Float64
    axis::A    
end

function FrameLocalVelocitySensor(;body::Body2D, position::SVector{2, Float64}, rot::Float64=0.0, axis::Union{Int64,Symbol})
    if axis == 1 || :x || axis == :X
        return FrameLocalVelocitySensor(body, position, rot, XAxis())
    elseif axis == 2 | :y || axis == :Y
        return FrameLocalVelocitySensor(body, position, rot, YAxis())
    else
        throw(ArgumentError("Axis must be either 1 (for X) or  2 (for Y) "))
    end
end


function measure(state::Vector{Float64}, dstate::Vector{Float64}, sys::MBSystem2D, sensor::FrameLocalVelocitySensor{XAxis})::Float64
    body = sensor.body;
    pos_dofs = get_body_position_dofs(sys, body);
    vel_dofs = get_body_velocity_dofs(sys, body);
    x_, y_, r_ = state[pos_dofs];
    vx_, vy_, vr_ = state[vel_dofs];
    s_x, s_y = sensor.position;
    s_r = sensor.rot;
    sin_α, cos_α = sincos(s_r);
    sin_θ, cos_θ = sincos(r_);

    return (vx_ * cos_θ + vy_ * sin_θ + vr_ * s_y) * cos_α + (-vx_ * sin_θ + vy_ * cos_θ + vr_ * s_x) * sin_α 
end
function measure(state::Vector{Float64}, dstate::Vector{Float64}, sys::MBSystem2D, sensor::FrameLocalVelocitySensor{YAxis})::Float64
    body = sensor.body;
    pos_dofs = get_body_position_dofs(sys, body);
    vel_dofs = get_body_velocity_dofs(sys, body);
    x_, y_, r_ = state[pos_dofs];
    vx_, vy_, vr_ = state[vel_dofs];
    s_x, s_y = sensor.position;
    s_r = sensor.rot;
    sin_α, cos_α = sincos(s_r);
    sin_θ, cos_θ = sincos(r_);

    return -(vx_ * cos_θ + vy_ * sin_θ - vr_ * s_y) * sin_α + (-vx_ * sin_θ + vy_ * cos_θ + vr_ * s_x) * cos_α
end

mutable struct FrameLocalAccelerationSensor{A<:AbstractAxis} <: AbstractSensor2D
    body::Body2D
    position::SVector{2, Float64}
    rot::Float64
    axis::A    
end

function FrameLocalAccelerationSensor(;body::Body2D, position::SVector{2, Float64}, rot::Float64=0.0, axis::Union{Int64, Symbol})
    if axis == 1 || :x || axis == :X
        return FrameLocalAccelerationSensor(body, position, rot, XAxis())
    elseif axis == 2 || :y || axis == :Y
        return FrameLocalAccelerationSensor(body, position, rot, YAxis())
    else
        throw(ArgumentError("Axis must be either 1 (for X) or  2 (for Y) "))
    end
end


function measure(state::Vector{Float64}, dstate::Vector{Float64}, sys::MBSystem2D, sensor::FrameLocalAccelerationSensor{XAxis})::Float64
    body = sensor.body
    pos_dofs = get_body_position_dofs(sys, body)
    vel_dofs = get_body_velocity_dofs(sys, body)
    
    θ = state[pos_dofs[3]]
    ω = state[vel_dofs[3]]
    
    ax_cm = dstate[vel_dofs[1]]
    ay_cm = dstate[vel_dofs[2]]
    α_ang = dstate[vel_dofs[3]]
    
    sx, sy = sensor.position[1], sensor.position[2]
    
    s_θ, c_θ = sincos(θ)
    a_body_x =  ax_cm * c_θ + ay_cm * s_θ - α_ang * sy - (ω^2) * sx
    a_body_y = -ax_cm * s_θ + ay_cm * c_θ + α_ang * sx - (ω^2) * sy
    
    s_α, c_α = sincos(sensor.rot)
    return a_body_x * c_α + a_body_y * s_α
end
function measure(state::Vector{Float64}, dstate::Vector{Float64}, sys::MBSystem2D, sensor::FrameLocalAccelerationSensor{YAxis})::Float64
    body = sensor.body
    pos_dofs = get_body_position_dofs(sys, body)
    vel_dofs = get_body_velocity_dofs(sys, body)
    
    θ = state[pos_dofs[3]]
    ω = state[vel_dofs[3]]
    
    ax_cm = dstate[vel_dofs[1]]
    ay_cm = dstate[vel_dofs[2]]
    α_ang = dstate[vel_dofs[3]]
    
    sx, sy = sensor.position[1], sensor.position[2]
    
    s_θ, c_θ = sincos(θ)
    a_body_x =  ax_cm * c_θ + ay_cm * s_θ - α_ang * sy - (ω*ω) * sx
    a_body_y = -ax_cm * s_θ + ay_cm * c_θ + α_ang * sx - (ω*ω) * sy
    
    s_α, c_α = sincos(sensor.rot)
    return -a_body_x * s_α + a_body_y * c_α
end
