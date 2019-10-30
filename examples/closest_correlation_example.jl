include("../../ApproximateCOSMO.jl/src/COSMO.jl")
using Main.COSMO, Random
include("../src/closest_correlation.jl"); include("../utils/plot_utils.jl")
using RandomCorrelationMatrices

n = 300
C = RandomCorrelationMatrices.randcormatrix(n, 0.1) + 0.1*Symmetric(randn(n, n))
all(eigvals(C) .>= 0) && @warn "The perturbed correlation matrix is still positive definite!"

settings = COSMO.Settings(verbose=true, verbose_timing=true,
    eps_abs=1e-7, eps_rel=1e-7, max_iter=5000,
    rho=.1, adaptive_rho=false,
    check_termination=40, check_infeasibility=1000)
    # kkt_solver=COSMO.with_options(COSMO.IndirectReducedKKTSolver, solver=:CG)
model = create_closest_correlation_model(C, settings; psdcone_type=COSMO.PsdConeTriangleLanczos)
COSMO.optimize!(model)
