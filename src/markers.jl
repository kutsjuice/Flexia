mutable struct Marker <: AbstractMarker2D
    name::String
    body_index::Int64
    point_on_body::AbstractVector{Float64}
    
    function first_marker(;body::Body2D, sys::MBSystem2D, dofs::AbstractVector{Float64})
        name = "bmarker"
        body_index = get_index(body::Body2D)
        point_on_body = get_upperight_point(sys::MBSystem2D, body::Body2D, dofs::AbstractVector{Float64})
    end

    function second_marker(;body::Body2D, sys::MBSystem2D, dofs::AbstractVector{Float64})
        name = "bmarker"
        body_index = get_index(body::Body2D)
        point_on_body = get_downleft_point(sys::MBSystem2D, body::Body2D, dofs::AbstractVector{Float64})
    end
end


