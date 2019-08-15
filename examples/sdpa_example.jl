include("../../COSMO_original/src/COSMO.jl")
using Main.COSMO, JLD2
include("../src/sdpa.jl"); include("../utils/plot_utils.jl")

# Loading directly the SDPA datafile from a jld2 file. To convert files from the "dat-s" format to .jld2 see the script
# ../data/sdplib/sdpLibImport.jl
filename = "../data/sdpa/gpp500-1.jld2"
@load filename m nblocks blocks c F minimum_value

settings = COSMO.Settings(verbose=true, verbose_timing=true,
    eps_abs=1e-3, eps_rel=1e-3, max_iter=20000,
    rho=.1, adaptive_rho=false,
    check_termination=1000, check_infeasibility=40)
    #kkt_solver=COSMO.with_options(COSMO.IndirectReducedKKTSolver, solver=:CG)

model = create_sdpa_model(c, F, blocks, settings; psdcone_type=COSMO.PsdConeTriangleLanczos)
res = COSMO.optimize!(model)

using ProxSDP, MosekTools
solve_sdpa_jump(c, F, blocks, solver=ProxSDP.Optimizer)# plot_lanczos(model)