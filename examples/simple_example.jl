using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff
using JSON
using LinearAlgebra

using StaticArrays

const g = 9.81

bd1 = Body2D(1, 1)
bd2 = Body2D(1, 1)

bd2.forces[2] = (x) -> -bd2.mass * g

jnt1 = FixedJoint(bd1)
jnt2 = HingeJoint(bd1, bd2)

set_position_on_second_body!(jnt2, SA[-1., 0])

sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)

add!(sys, jnt1)
add!(sys, jnt2)

if (!assemble!(sys))
    println("Assembling failed!")
end

func = sys.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)

bd1_x_ind, bd1_y_ind, _ = get_body_position_dofs(sys, bd1)
bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(sys, bd2)

initial = zeros(number_of_dofs(sys))

func(initial)
jacoby(initial)

mass = zeros(number_of_dofs(sys), number_of_dofs(sys));
for i in 1:last_body_dof(sys)
    mass[i, i] = 1
end

time_start = 0
time_end = 10
time_step = 1

time_span = 0:0.001:10

sol2 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
cros!(sol2, initial, mass, func, jacoby, step(time_span))
animate(sys, sol2, time_span, "simple_bar.mp4"; framerate = 60, limits = (-2,13, -2, 13))

function save_vector(vector::Vector, filename::String)
    try
        open(filename, "w") do io
            JSON.print(io, vector)
        end
        println("✓ Vector saved to $filename")
    catch e
        println("✗ Error saving: $e")
    end
end

save_vector(sol2[:, end], "sol22_simple.json")

function eigenvector_matrix(A::AbstractMatrix)
    # Validate input is square
    m, n = size(A)
    if m != n
        error("Input matrix must be square, got $(m)×$(n)")
    end
    
    # Compute eigen decomposition
    F = eigen(A)
    
    # Return eigenvector matrix (each column is an eigenvector)
    return F.vectors
end

function get_imaginary_part(matrix::AbstractMatrix)
    return imag.(matrix)
end


eigen_matrix = get_imaginary_part(eigenvector_matrix(jacoby(sol2[:,end])))

function save_matrix_to_json(matrix::AbstractMatrix, filename::String)
    data = Dict{String, Any}(
        "eltype" => string(eltype(matrix)),
        "size" => collect(size(matrix)),
        "data" => matrix
    )
    
    open(filename, "w") do io
        JSON.print(io, data, 2)
    end
    
    println("✓ Matrix saved to $filename")
    return true
end

function write_matrix_formatted(filename::String, matrix::AbstractMatrix)
    open(filename, "w") do file
        # Записываем открывающую скобку
        write(file, "[")
        
        nrows, ncols = size(matrix)
        
        for i in 1:nrows
            # Переходим на новую строку для всех строк кроме первой
            if i > 1
                write(file, " ")
            end
            
            # Записываем элементы строки
            for j in 1:ncols
                write(file, string(matrix[i, j]))
                if j < ncols
                    write(file, " ")
                end
            end
            
            # Записываем разделитель строк или закрывающую скобку
            if i < nrows
                write(file, "\n")
            else
                write(file, "]")
            end
        end
    end
end

save_matrix_to_json(eigen_matrix, "simple_eigen_matrix.json")
write_matrix_formatted("simple_eigen_matrix.txt", eigen_matrix)
s = number_of_dofs(sys)
println("$s")
write_matrix_formatted("simple_solve.txt", sol2)