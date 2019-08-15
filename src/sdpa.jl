using SparseArrays
using Main.COSMO


function create_sdpa_model(c::Vector{T}, F, blocks, settings=nothing; psdcone_type=COSMO.PsdConeTriangle)  where {T}
    m = length(c)
    constraints = COSMO.Constraint{T}[]
    for block_id in 1:length(blocks)
        block = blocks[block_id]
        if block > 0
            if psdcone_type != COSMO.PsdCone
                A = hcat([COSMO.extract_upper_triangle(F[i+1, block_id]) for i = 1:m]...)
            else
                A = hcat([reshape(F[i+1, block_id], length(F[i+1, block_id])) for i = 1:m]...)
            end

            if psdcone_type != COSMO.PsdCone
                b = -COSMO.extract_upper_triangle(Matrix(F[1, block_id]))
            else
                b = -Vector(reshape(F[1, block_id], length(F[1, block_id])))
            end
            push!(constraints, COSMO.Constraint(A, b, psdcone_type))
        else
            A = hcat([diag(F[i+1, block_id]) for i = 1:m]...)
            b = -Vector(diag(F[1, block_id]))
            push!(constraints, COSMO.Constraint(A, b, COSMO.Nonnegatives))
        end
    end
    P = SparseMatrixCSC(zero(T)*I, m, m) # Hessian of the objective

    model = COSMO.Model()
    if settings != nothing
        COSMO.assemble!(model, P, c, constraints, settings=settings)
    else
        COSMO.assemble!(model, P, c, constraints)
    end
    return model
end


using JuMP, SCS

function solve_sdpa_jump(c, F, blocks; solver=SCS.Optimizer, kwargs...)
    m = length(c)
    model = Model(with_optimizer(solver, log_verbose=true))
    @variable(model, x[1:m])
    for block_id in 1:length(blocks)
        block_dimension = abs(blocks[block_id])
        
        if block_dimension > 1 
            b = -Vector(reshape(F[1, block_id], length(F[1, block_id])))
            A = hcat([reshape(F[i+1, block_id], length(F[i+1, block_id])) for i = 1:m]...)
            @constraint(model, A*x + b in MOI.PositiveSemidefiniteConeSquare(block_dimension))
        else # Scalar case
            b = -Vector(diag(F[1, block_id]))
            A = hcat([diag(F[i+1, block_id]) for i = 1:m]...)
            @constraint(model, A*x + b >= 0)
        end
    end
    @objective(model, Min, dot(c, x))
    JuMP.optimize!(model)
    @show JuMP.termination_status(model)
    return value.(x), JuMP.objective_value(model)
end