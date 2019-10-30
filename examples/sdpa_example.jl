include("../../ApproximateCOSMO.jl/src/COSMO.jl")
using Main.COSMO, JLD2
include("../src/sdpa.jl")
using Printf
using DataFrames, CSV
using Glob
using Base.Filesystem
# include("../utils/plot_utils.jl")

@show LinearAlgebra.BLAS.openblas_get_config()

filepaths = glob("../data/sdpa/*.dat-s")

df = DataFrame(name=String[], time = Float64[], time_lobpcg = Float64[],
    time_prj = Float64[], time_prj_lobpcg = Float64[],
    iter = Int[], iter_lobpcg = Int[],
    objective = Float64[], objective_lobpcg = Float64[],
    optimal_objective = Float64[], status = [], status_lobpcg = []
)

filepaths = [filepaths[1]; filepaths]
# filepaths = filepaths[20:end]
for filepath in filepaths[1:end]
    # if basename(filepath)[1:end-6] == "maxG55" #&& !in(basename(filepath)[1:end-6], ["maxG60" "maxG55" "qpG51" "maxG32"])
    c, F, A, b, optimal_objective = load_sdpa_file(filepath)

    _, _, objective_lobpcg, status_lobpcg, solver_lobpcg = solve_sdpa_jump_dual(c, F, A, b, COSMO.Optimizer,
        psd_projector = COSMO.PsdConeTriangleLOBPCG,
        eps_abs = 1e-4, eps_rel = 1e-4,
        verbose = true,
        adaptive_rho_tolerance = 10.0,
        # max_iter = 5000000, check_termination = 100000,
        # adaptive_rho = false, rho=1.2
    )
    t_lobpcg = solver_lobpcg.times.solver_time
    t_proj_lobpcg = solver_lobpcg.times.proj_time
    iterations_lobpcg = solver_lobpcg.iterations

    _, _, objective, status, solver = solve_sdpa_jump_dual(c, F, A, b, COSMO.Optimizer,
        eps_abs = 1e-4, eps_rel = 1e-4,
        verbose = true,
        adaptive_rho_tolerance = 10.0,
        # max_iter = 5000000, check_termination = 100000,
        # adaptive_rho = false, rho=1.2
    )
    t = solver.times.solver_time
    t_proj = solver.times.proj_time
    iterations = solver.iterations
    #=
    status = NaN
    objective = NaN
    t = NaN
    t_proj = NaN
    iterations = -1
    =#

    push!(df, [basename(filepath)[1:end-6], t, t_lobpcg,
    t_proj, t_proj_lobpcg,
    iterations, iterations_lobpcg,
     objective, objective_lobpcg, optimal_objective,
     string(status), string(status_lobpcg)])
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
