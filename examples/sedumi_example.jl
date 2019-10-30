include("../../ApproximateCOSMO.jl/src/COSMO.jl")
using Main.COSMO
include("../src/sedumi.jl"); include("../utils/plot_utils.jl")
using SparseArrays, MAT

A, b, c, K = load_sedumi_file(filename)
x1, value1, status1 = solve_sedumi_jump(A, b, c, K, COSMO.Optimizer)
x2, value2, status2 = solve_sedumi_jump_dual(A, b, c, K, COSMO.Optimizer, max_iter = 500000, check_termination = 1000)
