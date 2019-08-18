include("../../COSMO_original/src/COSMO.jl")
using Main.COSMO, JLD2
include("../src/sdpa.jl"); include("../utils/plot_utils.jl")

# Load Objectives
using DataFrames, CSV
df = CSV.File("../data/sdpa/sdplib_info.csv") |> DataFrame
objectives = Dict(zip(df[1], df[4]))

file = "truss1"
@time c, F, A, b = load_sdpa_file(string("../data/sdpa/", file, ".dat-s")
objective = objectives[file]

solution, solution_objective = solve_sdpa_jump(c, F, A, b, COSMO.Optimizer)
solution, solution_objective = solve_sdpa_jump_dual(c, F, A, b, COSMO.Optimizer,
    max_iter = 500000, eps_abs = 1e-5, eps_rel = 1e-5, check_termination = 5000)