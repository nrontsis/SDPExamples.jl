using Plots

function plot_lanczos(model)
    plotly()
    i = 0;
    for set in model.p.C.sets
        if isa(set, COSMO.PsdConeTriangleLanczos)
            p1 = plot(hcat(set.subspace_dim_history,
                    set.n/2*ones(length(set.subspace_dim_history))),
                    title=string("Set #", i,":  Subspace sizes"),
                    label=["subspace size" "n/2"],
                    # line=:stem
                )
            p2 = plot(hcat(max.(set.residual_history, 1e-6),
                    max.(set.λ_rem_history, 1e-6)),
                    title=string("Set #", i,":  Residual history"),
                    yaxis=:log,
                    label=["total residual" "only from ignored subspace"],
                    # line=:stem
                )
            p3 = scatter(set.λ_rem_multiplications, 
                    title=string("Set #", i,":  ARPACK multiplications"),
                    markersize = 2,
                )
            plot(p1, p2, p3)
            gui()
        end
        i = i + 1
	end
end