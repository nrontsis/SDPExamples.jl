using JuMP, MAT

function load_sedumi_file(filename)
   	file = matopen(filename)
   	if exists(file, "A")
      		A = read(file, "A");
   	elseif exists(file, "At")
      		A = read(file, "At")';
   	else
      		@assert false "The MATLAB file must contain either A or At as a variable"
   	end
   	b = Vector(read(file, "b")[:]);
   	c = Vector(read(file, "c")[:]);
   	K = read(file, "K")
   	close(file)

	# Cast Int entries
   	for key in ["f", "l"]
      		if haskey(K, key)
         			K[key] = Int(K[key])
      		else
         			K[key] = 0
      		end
   	end
	# Cast Vector{Int} entries
   	for key in ["q", "r", "s"]
      		if haskey(K, key)
         			@assert isa(K[key], AbstractArray) || isa(K[key], Number)

         			if isa(K[key], AbstractArray)
            				K[key] = Vector{Int}(K[key][:])
         			else # isa(K[key], Number)
            				K[key] = [Int(K[key])]
         			end
         			K[key] = K[key][K[key] .> 0]
      		else
         			K[key] = Int[]
      		end
   	end

   	return A, b, c, K
end

function solve_sedumi_jump(A, b, c, K, solver; kwargs...)
	# Solve Problems in Sedumi form, i.e.
	# minimize    c'x
	# subject to  Ax = b
	#             x \in K
	# See https://github.com/sqlp/sedumi/blob/master/sedumi.m#L49-L92
	
   	model = Model(with_optimizer(solver; kwargs...))

	# Fisrt handle all variables except semidefinite matrices
   	constraint_types = [MOI.Reals;
		MOI.Nonnegatives;
		fill(MOI.SecondOrderCone, length(K["q"]));
		fill(MOI.RotatedSecondOrderCone, length(K["r"]))]
   	dimensions = [K["f"]; K["l"]; K["q"]; K["r"]]

   	x = Vector{VariableRef}(undef, 0)
   	for (dimension, constraint_type) in zip(dimensions, constraint_types)
      		x_ = @variable(model, [1:dimension])
      		if constraint_type != MOI.Reals # Some solvers don't accept MOI.Reals
         			@constraint(model, x_ in constraint_type(dimension))
      		end
      		append!(x, x_)
   	end

	# Now handle semidefinite matrices
   	for dimension in K["s"]
      		x_ = @variable(model, [1:dimension, 1:dimension], PSD)
      		append!(x, reshape(x_, length(x_)))
   	end

   	@constraint(model, A * x .== b)
   	@objective(model, Min, dot(c, x))

   	JuMP.optimize!(model)
    return value.(x), JuMP.objective_value(model), JuMP.termination_status(model)
end

function solve_sedumi_jump_dual(A, b, c, K, solver; kwargs...)
	# Solve the dual of problem described in Sedumi form, i.e.
	# maximize    b'y
	# subject to  c - A'y \in K*
	# See https://github.com/sqlp/sedumi/blob/master/sedumi.m#L49-L92

   	model = Model(with_optimizer(solver; kwargs...))

	# Fisrt handle all variables except semidefinite matrices
   	constraint_types = [MOI.Zeros;
		MOI.Nonnegatives;
		fill(MOI.SecondOrderCone, length(K["q"]));
		fill(MOI.RotatedSecondOrderCone, length(K["r"]))]
   	dimensions = [K["f"]; K["l"]; K["q"]; K["r"]]

   	@variable(model, y[1:length(b)])
   	start_idx = 1
   	for (dimension, constraint_type) in zip(dimensions, constraint_types)
      		indices = start_idx:start_idx + dimension - 1
      		dimension > 0 && @constraint(model, c[indices] - A[:, indices]' * y in constraint_type(dimension))
      		start_idx += dimension
   	end
   	for dimension in K["s"]
      		indices = Vector(start_idx:start_idx + dimension^2 - 1)[upper_triangular_indices(dimension)]
      		@constraint(model, c[indices] - A[:, indices]' * y in MOI.PositiveSemidefiniteConeTriangle(dimension))
      		start_idx += dimension^2
   	end

   	@objective(model, Max, dot(b, y))

   	JuMP.optimize!(model)
    return value.(y), JuMP.objective_value(model), JuMP.termination_status(model)
end

function upper_triangular_indices(n)
   	indices = Vector{Int}(undef, Int(n * (n + 1) / 2))
   	index = 0
   	counter = 0
   	for j = 1:n, i = 1:n
      		index += 1
      		if i <= j
         			counter += 1
         			indices[counter] = index
      		end
   	end
   	return indices
end