

using Main.COSMO
using Test, LinearAlgebra, SparseArrays, Random

function create_diagonal_extractor(Tv, n::Ti) where {Ti}
    # Creates a matrix A that slices out the diagonal entries Xii of
    # a upper-triangular-vectorized square matrix x=vec(X) of dimension n
    A = spzeros(Tv, Ti, n, Ti(n*(n + 1)/2))
    idx = 0
    for i = 1:n
        idx += i
        A[i, idx] = one(Tv)
    end
    return A
end

function create_closest_correlation_model(C::AbstractMatrix{Tv}, settings=nothing;
    psdcone_type=COSMO.PsdConeTriangle) where Tv
    # Closest Correlation Matrix problem
    # min   1/2 ||X-C||^2
    # s.t.    Xii = 1
    #         X ⪴ 0

    @assert psdcone_type in [COSMO.PsdConeTriangle COSMO.PsdConeTriangleLanczos]
    n = size(C, 1)
    @assert size(C, 2) == size(C, 1)

    # Contraint #1: Diagonal equal to one
    A1 = create_diagonal_extractor(Tv, n)
    b1 = -ones(Tv, n)
    constraints = [COSMO.Constraint(A1, b1, COSMO.ZeroSet)]

    # Contraint #2: decision variable is a PSD
    n2 =  Int(n*(n + 1)/2)
    A2 = SparseMatrixCSC(one(Tv)*I, n2, n2)
    b2 = zeros(Tv, n2)
    push!(constraints, COSMO.Constraint(A2, b2, psdcone_type))

    P = SparseMatrixCSC(Tv(1/2)*I, n2, n2)
    q = -COSMO.extract_upper_triangle(C)

    model = COSMO.Model()
    if settings != nothing
        COSMO.assemble!(model, P, q, constraints, settings=settings)
    else
        COSMO.assemble!(model, P, q, constraints)
    end

    return model
end
