include("../../COSMO_original/src/COSMO.jl")
using Main.COSMO, JLD2
include("../src/sdpa.jl"); include("../utils/plot_utils.jl")
using Printf
using DataFrames, CSV
using Glob
using Base.Filesystem

filepaths = glob("../data/sdpa/*.dat-s")

df = DataFrame(name=String[], time = Float64[], time_lanczos = Float64[],
    objective = Float64[], objective_lanczos = Float64[],
    optimal_objective = Float64[], status = [], status_lanczos = []
)

# filepath = "../data/sdpa/equalG51.dat-s"
for filepath in filepaths[1:end]
    c, F, A, b, optimal_objective = load_sdpa_file(filepath)
    objective_lanczos = nothing; t_lanczos = Inf; status_lanczos = nothing
    for i = 1:2
        _, _, objective_lanczos, status_lanczos, t_lanczos = solve_sdpa_jump_dual(c, F, A, b, COSMO.Optimizer,
            lanczos = true,
            eps_abs = 1e-4, eps_rel = 1e-4,
            # max_iter = 5000000, check_termination = 100000,
            # adaptive_rho = false, rho=1.2
        )
    end
    objective = nothing; t = Inf; status = nothing
    for i = 1:2
        _, _, objective, status, t = solve_sdpa_jump_dual(c, F, A, b, COSMO.Optimizer,
            lanczos = false,
            eps_abs = 1e-4, eps_rel = 1e-4,
            # max_iter = 5000000, check_termination = 100000,
            # adaptive_rho = false, rho=1.2
        )
    end
    push!(df, [basename(filepath)[1:end-5], t, t_lanczos,
     objective, objective_lanczos, optimal_objective,
     string(status), string(status_lanczos)])
    df |> CSV.write(string("results.csv"))
end

#=
using MosekTools
solution, solution_objective = solve_sdpa_jump_dual(c, F, A, b, Mosek.Optimizer)

using ProxSDP
for i = 1:2
    solution, solution_objective = solve_sdpa_jump_dual(c, F, A, b, ProxSDP.Optimizer, log_verbose=true,
        tol_primal = 1e-4, tol_dual = 1e-4)
end
=#