using LinearAlgebra, SparseArrays
using JuMP
using CSV, DataFrames
using Base.Filesystem
using Printf
using JLD2

function load_sdpa_file(filepath)
    info = CSV.File(joinpath(dirname(filepath), "sdplib_info.csv"))  |> DataFrame
    objectives = Dict(zip(info[1], info[4]))
    filename =  join(split(basename(filepath), ".")[1:end-1], ".")
    objective = haskey(objectives, filename) ? objectives[filename] : NaN
    @printf("Problem: %s with optimal objective: %.4e\n", filename, objective)
    # Load Problem
    jld_filepath = join([split(filepath, ".")[1:end-1]; "jld2"], ".")
    if isfile(jld_filepath) # Attempt to load saved jld2 file
        @load jld_filepath c F A b
    else # Load directly from dat-s file
        print("Loading SDPA file. This might take some time... ")
        t = @elapsed c, F, A, b = _load_sdpa_file(string(jld_filepath[1:end-5], ".dat-s"))
        println(string("Done in ", t, " seconds!"))
        @save jld_filepath c F A b
    end
    return c, F, A, b, objective
end

function _load_sdpa_file(filename)
    #=
    Load problems of the SDPA format from "dat-s" files, i.e. problems of the form
    
    minimize c'x
    subject to  F1*x1 + ... + Fm*xm - F0 ⪰ 0

    with the real m-dimensional vector x as the decision variable.
    See https://freesoft.dev/program/129366650#sdpa-sparse-format for a detailed description.

    Note that the coefficient matrices Fi have common block diagonal structure. Furthermore, some of 
    these blocks are themselves diagonal. Thus we can transform the original problem into

    minimize c'x
    subject to  F1,1*x1 + ... + Fm,1*xm - F0,1 ⪰ 0
                ...
                F1,l*x1 + ... + Fm,l*xm - F0,l ⪰ 0
                A*x ≥ b

    where Fi,j denotes the j-th diagonal block of the i-th coefficient matrix and the linear inequality constraints
    arised from subblocks of the coefficient matrices that have diagonal structure.

    The function returns
    - c::Vector the linear cost
    - F::Array{SparseMatrixCSC,2} F[i, j] is the j-th block of the (i-1)-th coefficient matrix
    - A::SparseMatrixCSC, b::Vector: for the linear inequality constraints
    - opt_val::Float64 the optimal value of the problem. Only supported from problems in the SDPLIB collection.
        These are extracted from a README file of the SDPLIB suite.
        Contains NaN if README is not present or if the extraction of the optimal value failed.
    =#
    file = open(filename);
    lines = readlines(file)
    close(file)

    # Discard all lines that contain comments
    indices = findall(line->length(line) == 0 || (line[1] != '*' && line[1] != '"'), lines)
    lines = lines[indices]

    m = parse(Int, strip(lines[1])) # length of decision variable
    blocks_sizes_string = split(replace(strip(lines[3], ['{','}','(',')']), "+" => ""))
    block_sizes = [parse(Int, str) for str in blocks_sizes_string]
    @assert length(block_sizes) == parse(Int, strip(lines[2]))

    # try first to split by whitespace as delimiter (also remove certain trailing, ending characters)
    c_strings_vector = split(replace(strip(lines[4], ['{','}','(',')']), "+" => ""))
    # otherwise try comma as delimiter
    if length(c_strings_vector) != m
        c_strings_vector = split(replace(strip(lines[4], ['{','}','(',')']), "+" => ""), ",")
    end
    c = [parse(Float64, str) for str in c_strings_vector]
    
    @assert length(c) == m
    indices = Vector(1:length(block_sizes));
    indices_map = [indices[block_sizes .> 1]; indices[block_sizes .<= 1]]
    diag_sizes = abs.(block_sizes[block_sizes .<= 1])
    block_sizes = block_sizes[block_sizes .> 1]

    F = [spzeros(block_sizes[j], abs(block_sizes[j])) for i = 1:m + 1, j = 1:length(block_sizes)]
    A = spzeros(sum(diag_sizes), m)
    b = zeros(size(A, 1))
    for line in lines[5:end]
        parsed_line = [parse(Float64, str) for str in split(line)]
        matrix_number = Int(parsed_line[1]) + 1
        block_number = findfirst(indices_map .== Int(parsed_line[2]))
        
        i = Int(parsed_line[3])
        j = Int(parsed_line[4])
        entry = parsed_line[5]

        if block_number <= length(block_sizes)
            F[matrix_number, block_number][i, j] = entry
            F[matrix_number, block_number][j, i] = entry
        else
            # Diagonal block - convert to linear inequalities
            @assert i == j
            idx = i + sum(diag_sizes[1:block_number - length(block_sizes) - 1])

            if matrix_number == 1
                b[idx] = entry
            else
                A[idx, matrix_number - 1] = entry
            end
        end 
    end

    return c, F, A, b
end

function extract_upper_triangle(A::SparseMatrixCSC{Tv,Ti}, scaling_factor::Tv = one(Tv)) where {Tv,Ti}
   	result_nnz = Ti(nnz(A) / 2 + nnz(diag(A)) / 2)
   	nzind = zeros(Ti, result_nnz)
   	nzval = zeros(Tv, result_nnz)
   	n = size(A, 1)
    counter = 0
   	for j in 1:n, idx in A.colptr[j]:A.colptr[j + 1] - 1
        i = A.rowval[idx]
      		if i <= j
            counter += 1
            nzind[counter] = Ti(j * (j - 1) / 2 + i)
            if i == j
                nzval[counter] = A.nzval[idx]
            else
                nzval[counter] = scaling_factor * A.nzval[idx]
         			end
      		end
    end
   	return SparseVector{Tv,Ti}(Ti(n * (n + 1) / 2), nzind, nzval)
end

function solve_sdpa_jump(c, F, A, b, solver; kwargs...)
    m = length(c)
    model = Model(with_optimizer(solver; kwargs...))
    @variable(model, x[1:m])
    @constraint(model, A * x .>= b)
    for block_id in 1:size(F, 2)
        A_ = hcat([extract_upper_triangle(F[i + 1, block_id]) for i = 1:m]...)
        b_ = Vector(extract_upper_triangle(F[1, block_id]))
        block_size = size(F[1, block_id], 1)
        @constraint(model, A_ * x - b_ in MOI.PositiveSemidefiniteConeTriangle(block_size))
    end
    @objective(model, Min, dot(c, x))
    JuMP.optimize!(model)
    raw_solver = MOI.get(model, MOI.RawSolver())
    return value.(x), JuMP.objective_value(model), JuMP.termination_status(model), raw_solver
end

function solve_sdpa_jump_dual(c, F, A, b, solver; kwargs...)
    model = Model(with_optimizer(solver; kwargs...))

    @variable(model, y[1:size(A, 1)] >= 0)
    num_blocks = size(F, 2)
    block_sizes = [size(F[1, block_id], 1) for block_id in 1:num_blocks]
    Y = [@variable(model, [1:block_size, 1:block_size], PSD) for block_size in block_sizes]

    constraint_sum = A' * y
    objective_sum = dot(y, b)
    for block_id in 1:num_blocks
        # Normally, we should also avoid the matrix dot product below, but in most solvers this doesn't matter
        add_to_expression!.(objective_sum, dot(F[1, block_id], Y[block_id]))
        for i = 1:length(c)
            f = SparseVector(F[i + 1, block_id][:])
            add_to_expression!(constraint_sum[i], dot(f.nzval, Y[block_id][f.nzind]))
        end
    end
    @constraint(model, constraint_sum .== c)
    @objective(model, Max, objective_sum)
    
    JuMP.optimize!(model)
    raw_solver = MOI.get(model, MOI.RawSolver())

    return [value.(Y[i]) for i = 1:length(block_sizes)], value.(y),
        JuMP.objective_value(model),
        JuMP.termination_status(model),
        raw_solver
end

#=
function solve_sdpa_jump_b(c, F, A, b; solver=SCS.Optimizer, kwargs...)
    # A version of solve_sdpa_jump_b that does not explicitly construct A_, b_
    # It is significantly slower than the original function, solve_sdpa_jump.
    m = length(c)
    model = Model(with_optimizer(solver; kwargs...))
    @variable(model, x[1:m])
    @constraint(model, A*x .>= b)
    for block_id in 1:size(F, 2)
        block_size = size(F[1, block_id], 1)
        X = @variable(model, [1:block_size, 1:block_size], PSD) .- F[1, block_id]
        for i = 1:m
            add_to_expression!.(X, x[i].*F[i + 1, block_id])
        end
        @constraint(model, X .== 0)
    end
    @objective(model, Min, dot(c, x))
    JuMP.optimize!(model)
    return value.(x), JuMP.objective_value(model), JuMP.termination_status(model)
end
=#