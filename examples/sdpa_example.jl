include("../../COSMO_original/src/COSMO.jl")
using Main.COSMO, JLD2
include("../src/sdpa.jl")
using Printf
using DataFrames, CSV
using Glob
using Base.Filesystem
# include("../utils/plot_utils.jl")

filepaths = glob("../data/sdpa/*.dat-s")

df = DataFrame(name=String[], time = Float64[], time_lanczos = Float64[],
    time_prj = Float64[], time_prj_lanczos = Float64[],
    iter = Int[], iter_lanczos = Int[],
    objective = Float64[], objective_lanczos = Float64[],
    optimal_objective = Float64[], status = [], status_lanczos = []
)

# filepath = "../data/sdpa/equalG51.dat-s"
# filepaths = filepaths[20:end]
filepahts = [filepaths[1]; filepaths]
for filepath in filepaths[1:end]
    # if basename(filepath)[1:end-6] == "thetaG51" || !in(basename(filepath)[1:end-6], ["maxG60" "maxG55" "qpG51" "maxG32"])
    c, F, A, b, optimal_objective = load_sdpa_file(filepath)

    _, _, objective_lanczos, status_lanczos, solver_lanczos = solve_sdpa_jump_dual(c, F, A, b, COSMO.Optimizer,
        lanczos = true,
        eps_abs = 1e-4, eps_rel = 1e-4,
        verbose=false
        # max_iter = 5000000, check_termination = 100000,
        # adaptive_rho = false, rho=1.2
    )
    t_lanczos = solver_lanczos.times.solver_time
    t_proj_lanczos = solver_lanczos.times.proj_time
    iterations_lanczos = solver_lanczos.iterations

    _, _, objective, status, solver = solve_sdpa_jump_dual(c, F, A, b, COSMO.Optimizer,
        lanczos = false,
        eps_abs = 1e-4, eps_rel = 1e-4,
        verbose=false
        # max_iter = 5000000, check_termination = 100000,
        # adaptive_rho = false, rho=1.2
    )
    t = solver.times.solver_time
    t_proj = solver.times.proj_time
    iterations = solver.iterations

    push!(df, [basename(filepath)[1:end-5], t, t_lanczos,
    t_proj, t_proj_lanczos,
    iterations, iterations_lanczos,
     objective, objective_lanczos, optimal_objective,
     string(status), string(status_lanczos)])
    df |> CSV.write(string("results.csv"))
    flush(stdout)
    # end
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