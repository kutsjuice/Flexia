
function measure(state::Vector{Float64}, sys::MBSystem2D, sensor::AbstractSensor2D)::Float64
    return NaN64
end


mutable struct FramePositionSensor <: AbstractSensor2D
    body::Body2D
    position::SVector{2, Float64}
    rot::Float64
    crd::Int64
end
function FramePositionSensor(;body::Body2D, position::SVector{2, Float64}, rot::Number, crd::Int64)
    if crd < 1 || crd > 2
        throw("CRD is out of bounds")
    end
    return FramePositionSensor(body, position, rot, crd)
end
function FramePositionSensor(;body::Body2D, position::SVector{2, Float64}, rot::Number, crd::Symbol)
    crd_id = 0
    if crd == :xcord
        crd_id = 1
    elseif crd == :ycord
        crd_id = 2
    elseif crd == :rcord
        crd_id = 3
    end
    return FramePositionSensor(body, position, rot, crd_id)
end

function measure(state::Vector{Float64}, sys::MBSystem2D, sensor::FramePositionSensor)::Float64
    body = sensor.body;
    dofs = get_body_position_dofs(sys, body);
    x_, y_, r_ = state[dofs];
    if sensor.crd == 1
        return first(SA_F64[cos(r_) -sin(r_)] * sensor.position) + x_;
    else
        return first(SA_F64[sin(r_) cos(r_)] * sensor.position) + y_;
    end
end

function add!(sys::MBSystem2D, sensor::FramePositionSensor)
    push!(sys.sensors, sensor)
    return nothing
end

mutable struct FrameVelocitySensor <: AbstractSensor2D
end

mutable struct FrameAccelerationSensor <: AbstractSensor2D
end