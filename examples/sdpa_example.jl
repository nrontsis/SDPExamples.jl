include("../../COSMO_original/src/COSMO.jl")
using Main.COSMO, JLD2
include("../src/sdpa.jl"); include("../utils/plot_utils.jl")

@time c, F, A, b, minimum_value = load_sdpa_file("../data/sdpa/truss1.dat-s")
solution, solution_objective = solve_sdpa_jump(c, F, A, b, COSMO.Optimizer)
solution, solution_objective = solve_sdpa_jump_dual(c, F, A, b, COSMO.Optimizer,
    max_iter = 500000, eps_abs = 1e-5, eps_rel = 1e-5, check_termination = 5000)