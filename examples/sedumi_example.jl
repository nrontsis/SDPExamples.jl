include("../../COSMO_original/src/COSMO.jl")
using Main.COSMO
include("../src/sedumi.jl"); include("../utils/plot_utils.jl")
using SparseArrays, MAT


filename = "../data/sedumi/ks_D2L30N240.mat"
f = matopen(filename)
A = read(f, "A"); b = read(f, "b"); c = read(f, "c"); K = read(f, "K")
close(f)

settings = COSMO.Settings(verbose=true, verbose_timing=true,
    eps_abs=1e-4, eps_rel=1e-4, max_iter=50000,
    rho=.1, adaptive_rho=false,
    check_termination=1000, check_infeasibility=1000)
    # kkt_solver=COSMO.with_options(COSMO.IndirectReducedKKTSolver, solver=:CG)
    
model = create_sedumi_model(A, b, c, K, settings, psdcone_type=COSMO.PsdConeTriangle)
res = COSMO.optimize!(model)