using Main.COSMO
using SparseArrays

function create_sedumi_model(A, b, c, K, settings=nothing; psdcone_type=COSMO.PsdConeTriangle)
	@assert psdcone_type in [COSMO.PsdConeTriangle COSMO.PsdConeTriangleLanczos]

	# Cast types
    K["f"] = Int(K["f"]) # Free variables
    K["l"] = Int(K["l"]) # Nonnegative variables
    K["q"] = Vector{Int}(K["q"][:]) # Lorentz cones
    K["r"] = Vector{Int}(K["r"][:]) # Rotated Lorentz cones
    K["s"] = Vector{Int}(K["s"][:]) # PSD Cones

	m = K["f"] + K["l"] + sum(K["q"]) + sum(K["r"]) # Number of all variables except PSD cones
	if psdcone_type != COSMO.PsdCone 
		T = SparseMatrixCSC(I, m, m)
		for size_psd in (K["s"])
			T_ = get_triangulization_map(size_psd)
			T = [T 									spzeros(size(T, 1), size(T_, 2));
				spzeros(size(T_, 1), size(T, 2)) 	T_]
		end
	else
		T = I
	end
	AT = A*T'
	constraints = [COSMO.Constraint(AT, -b, COSMO.ZeroSet),]
	n = size(AT, 2)

    idx = 1
    if K["f"] > 0 # Free variables
        idx += K["f"]
    end
    if K["l"] > 0 # Nonnegative variables
        A_ = variable_selector(idx, K["l"], n)
        idx += size(A_, 1)
        push!(constraints, COSMO.Constraint(A_, zeros(size(A_, 1)), COSMO.Nonnegatives))
    end
	@assert length(K["q"]) == 0 && length(K["r"]) == 0 # No support for Lorentz or Rotated Lonentz cones
	for dim in K["s"] # Psd Variables
		if psdcone_type != COSMO.PsdCone 
			dim_unrolled = Int(dim*(dim + 1)/2)
		else
			dim_unrolled = n^2
		end
        A_ = variable_selector(idx, dim_unrolled, n)
        idx += size(A_, 1)
        push!(constraints, COSMO.Constraint(A_, -zeros(size(A_, 1)), psdcone_type))
    end

	model = COSMO.Model()
	P = SparseMatrixCSC(0.0*I, n, n) # Hessian of the objective
	# Linear objective
	if psdcone_type != COSMO.PsdCone 
		q = T*c
	else
		q = c
	end
	if settings != nothing
		COSMO.assemble!(model, P, q, constraints, settings=settings)
	else
		COSMO.assemble!(model, P, q, constraints)
	end

    return model
end

function variable_selector(start, dim, n)
	# Returns a SparseMatrixCSC that performs:
	# A*x = x[start:start+dim-1] where x is an n-dimensional vector
	@assert dim + start - 1 <= n "Variable selector our of bounds"
	A = spzeros(dim, n)
	@inbounds for i = 1:dim
		A[i, i + start - 1] = 1.0
	end
	return A
end

function get_upper_triangular_index(i, j, n)
	if i > j
		i, j = j, i
	end
	return Int((j - 1)*j/2 + i)
end

function get_triangulization_map(n)
	rowval = zeros(Int, n^2)
	nzval = zeros(Float64, n^2)
	colnum = n^2
	idx = 1
	for j = 1:n, i = 1:n
		rowval[idx] = get_upper_triangular_index(i, j, n)
		if i == j 
			nzval[idx] = 1.0
		else
			nzval[idx] = 1/sqrt(2)
		end

		idx += 1
	end
	colptr = 1 .+ [0; cumsum(ones(Int, colnum))]
	A = SparseMatrixCSC(Int(n*(n + 1)/2), n^2, colptr, rowval, nzval)
	return A
end


using JuMP, MosekTools

function solve_jump(A, b, c, K, solver=Mosek.Optimizer)
    n = length(c)
    K["f"] = Int(K["f"])
    K["r"] = Vector{Int}(K["r"][:])
    K["l"] = Int(K["l"])
    K["s"] = Vector{Int}(K["s"][:])
    K["q"] = Vector{Int}(K["q"][:])
    
    model = Model(with_optimizer(solver))

	# ToDo: FixMe
    @assert length(K["s"]) == 2  # This code is taylored to exactly two psd variables.
    @variable(model, x[1:K["f"]])
    @variable(model, X[1:K["s"][1], 1:K["s"][1]], PSD)
    @variable(model, Y[1:K["s"][2], 1:K["s"][2]], PSD)
    @constraint(model, A*[x; X[:]; Y[:]] .== b)

    @objective(model, Min, dot([x; X[:]; Y[:]], c))

    JuMP.optimize!(model)
    JuMP.objective_value(model)
    X = JuMP.value.(X)
    Y = JuMP.value.(Y)
    x = JuMP.value.(x)
    return X, Y, x, A, b, c, K
end