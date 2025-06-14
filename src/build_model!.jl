"""
    build_model!(mc; max_iter=50, constr_viol_tol=1e-8, bound_push=1e-15)

Builds the GTAP model in the provided model container (`mc`) which has an empty Ipopt model. This function is not typically called by the user; instead, it is invoked by `generate_initial_model()`. However, if the user wishes to change the model, he or she can extend this function.

### Args:
- `mc`: a model container

### Optional args:

- `max_iter`: maximum number of iterations
- `constr_viol_tol`: accuracy for constraint satisfaction
- `bound_push`: mandatory move of the starting values from constraint bounds


### Retruns:

Nothing but it modifies the `mc` object, adding on the existing `model` element
        
### Example:

`julia
build_model!(mc)
`
"""
function build_model!(mc; max_iter=50, constr_viol_tol=1e-8, bound_push=1e-15, calibration=true)
    # Structural parameters (some CES/CET options are not happening)
    δ_evfp = .!isnan.(mc.data["α_qfe"]) .&& mc.data["α_qfe"] .!= 0
    δ_maks = .!isnan.(mc.data["α_qca"]) .&& mc.data["α_qca"] .!= 0
    δ_vtwr = .!isnan.(mc.data["α_qtmfsd"]) .&& mc.data["α_qtmfsd"] .!= 0
    δ_vtwr_sum = NamedArray(.!mapslices(all, isnan.(mc.data["α_qtmfsd"]) .|| mc.data["α_qtmfsd"] .== 0, dims=1)[1, :, :, :], names(mc.data["vtwr"])[[2, 3, 4]])
    δ_qxs = .!isnan.(mc.data["α_qxs"]) .&& mc.data["α_qxs"] .!= 0
    δ_qga = .!isnan.(mc.data["α_qga"]) .&& mc.data["α_qga"] .!= 0
    δ_qia = .!isnan.(mc.data["α_qia"]) .&& mc.data["α_qia"] .!= 0
    δ_qfa = .!isnan.(mc.data["α_qfa"]) .&& mc.data["α_qfa"] .!= 0
    δ_act = .!isnan.(mc.data["γ_qca"]) .&& mc.data["γ_qca"] .!= 0


    # Read  sets
    (; reg, comm, marg, acts, endw, endwc, endws, endwm, endwms, endwf) = NamedTuple(Dict(Symbol(k) => mc.sets[k] for k ∈ keys(mc.sets)))

    # Read hard parameters
    (; esubt, esubc, esubva, esubd, etraq, esubq, subpar, incpar, etrae, esubg, esubm, esubs, rordelta, rorflex) = NamedTuple(Dict(Symbol(k) => mc.parameters[k] for k ∈ keys(mc.parameters)))

    # All variables used in the model
    @variables(mc.model,
        begin
            # Population
            pop[reg]

            # Firms top nest
            qint[acts, reg]
            pint[acts, reg]
            qva[acts, reg]
            pva[acts, reg]
            qo[acts, reg]
            po[acts, reg]

            # Firms second nest
            qfa[comm, acts, reg]
            pfa[comm, acts, reg]
            qfe[endw, acts, reg]
            pfe[endw, acts, reg]
            tfe[endw, acts, reg]
            qfd[comm, acts, reg]
            pfd[comm, acts, reg]
            qfm[comm, acts, reg]
            pfm[comm, acts, reg]

            # # # Firm distribution
            qca[comm, acts, reg]
            pca[comm, acts, reg]
            ps[comm, acts, reg]
            qc[comm, reg]
            pds[comm, reg]
            to[comm, acts, reg]

            # Endowments
            peb[endw, acts, reg]
            qes[endw, acts, reg]
            pfactor[reg]

            # Income
            fincome[reg]
            y[reg]
            p[reg]
            u[reg]
            ug[reg]

            # Private consumption        
            yp[reg]
            up[reg]
            uelas[reg]
            uepriv[reg]
            ppa[comm, reg]
            qpa[comm, reg]
            ppd[comm, reg]
            qpd[comm, reg]
            ppm[comm, reg]
            qpm[comm, reg]
            ppriv[reg]

            # Government consumption
            yg[reg]
            pgov[reg]
            pga[comm, reg]
            qga[comm, reg]
            pgd[comm, reg]
            qgd[comm, reg]
            pgm[comm, reg]
            qgm[comm, reg]

            # Saving
            qsave[reg]
            psave[reg]

            # Investment consumption
            pia[comm, reg]
            qia[comm, reg]
            pid[comm, reg]
            qid[comm, reg]
            pim[comm, reg]
            qim[comm, reg]
            qinv[reg]
            pinv[reg]

            # Trade - exports
            qms[comm, reg]
            qxs[comm, reg, reg]
            pmds[comm, reg, reg]
            pms[comm, reg]

            # Trade - margins
            ptrans[comm, reg, reg]
            qtmfsd[marg, comm, reg, reg]
            pt[marg]
            qtm[marg]
            qst[marg, reg]

            # Trade - imports/exports
            txs[comm, reg, reg]
            tx[comm, reg]
            pfob[comm, reg, reg]
            pcif[comm, reg, reg]
            tms[comm, reg, reg]
            tm[comm, reg]

            # Domestic market clearing
            qds[comm, reg]

            # Taxes
            tfd[comm, acts, reg]
            tfm[comm, acts, reg]
            tpd[comm, reg]
            tpm[comm, reg]
            tgd[comm, reg]
            tgm[comm, reg]
            tid[comm, reg]
            tim[comm, reg]

            #Factor Market
            pes[endw, acts, reg]
            pe[endwms, reg]
            qe[endwms, reg]
            qesf[endwf, acts, reg]
            tinc[endw, acts, reg]

            # Global Investment
            globalcgds
            pcgdswld
            walras_sup
            walras_dem
            pfactwld

            # Capital stocks
            kb[reg]
            ke[reg]

            rorg
            rore[reg]
            rorc[reg]
            rental[reg]

            # Soft parameters
            α_qintva[["int", "va"], acts, reg]
            γ_qintva[acts, reg]
            α_qfa[comm, acts, reg]
            γ_qfa[acts, reg]
            α_qfe[endw, acts, reg]
            γ_qfe[acts, reg]
            α_qfdqfm[["dom", "imp"], comm, acts, reg]
            γ_qfdqfm[comm, acts, reg]
            α_qca[comm, acts, reg]
            γ_qca[acts, reg]
            α_pca[comm, acts, reg]
            γ_pca[comm, reg]
            σyp[reg]
            σyg[reg]
            α_qga[comm, reg]
            γ_qga[reg]
            α_qia[comm, reg]
            γ_qia[reg]
            α_qpdqpm[["dom", "imp"], comm, reg]
            γ_qpdqpm[comm, reg]
            α_qgdqgm[["dom", "imp"], comm, reg]
            γ_qgdqgm[comm, reg]
            β_qpa[comm, reg]
            α_qidqim[["dom", "imp"], comm, reg]
            γ_qidqim[comm, reg]
            α_qxs[comm, reg, reg]
            γ_qxs[comm, reg]
            α_qtmfsd[marg, comm, reg, reg]
            α_qst[marg, reg]
            γ_qst[marg]
            α_qes2[endws, acts, reg]
            γ_qes2[endws, reg]
            α_qinv[reg]

            δ[reg]
            ρ[reg]

        end
    )

    # GTAP has a weird thing where a parameter changes the model structure; we follow it here
    if rordelta==1
        @variables(mc.model,
        begin
        end
        )
    end

    if calibration
        @variables(mc.model,
            begin
                # Values
                vdpp[comm, reg]

                # Parameter constraints
                ϵ_qintva[acts, reg]
                ϵ_qfa[acts, reg]
                ϵ_qfe[acts, reg]
                ϵ_qfdqfm[comm, acts, reg]
                ϵ_qga[reg]
                ϵ_qia[reg]
                ϵ_qgdqgm[comm, reg]
                ϵ_qxs[comm, reg]
                ϵ_qes2[endws, reg]
                ϵ_qpdqpm[comm, reg]
                ϵ_qidqim[comm, reg]
                ϵ_qinv
                σsave[reg]

                # Shares (helpers to calibrate)
                σ_vp[comm, reg]
                σ_vdp[comm, reg]
                σ_vg[comm, reg]
                σ_vdg[comm, reg]
                σ_vi[comm, reg]
                σ_vdi[comm, reg]
                σ_vf[comm, acts, reg]
                σ_vdf[comm, acts, reg]
                σ_vff[endw, acts, reg]
                σ_vif[acts, reg]
                σ_vtwr[marg, comm, reg, reg]
                σ_qxs[comm, reg, reg]
                σ_qinv[reg]
                σ_ρ[reg]
            end
        )
    end

    # Remove structurally missing variables

    delete.(mc.model, Array(γ_pca)[esubq.==0])

    # Remove missing activities
    delete.(mc.model, Array(qint[acts, reg])[.!δ_act])
    delete.(mc.model, Array(pint[acts, reg])[.!δ_act])
    delete.(mc.model, Array(qva[acts, reg])[.!δ_act])
    delete.(mc.model, Array(pva[acts, reg])[.!δ_act])
    delete.(mc.model, Array(qo[acts, reg])[.!δ_act])
    delete.(mc.model, Array(po[acts, reg])[.!δ_act])
    delete.(mc.model, Array(α_qintva["int",acts, reg])[.!δ_act])
    delete.(mc.model, Array(α_qintva["va",acts, reg])[.!δ_act])
    delete.(mc.model, Array(γ_qfa[acts, reg])[.!δ_act])
    delete.(mc.model, Array(γ_qfe[acts, reg])[.!δ_act])
    delete.(mc.model, Array(γ_qca[acts, reg])[.!δ_act])
    delete.(mc.model, Array(γ_qintva[acts, reg])[.!δ_act])
        

    if calibration
        delete.(mc.model, Array(ϵ_qintva)[.!δ_act])
        delete.(mc.model, Array(ϵ_qfa)[.!δ_act])
        delete.(mc.model, Array(σ_vif)[.!δ_act])
        delete.(mc.model, Array(ϵ_qfe)[.!δ_act])
    end

    # Remove qfe, qes, pfe, peb, α_qes2, α_qfe, cevos, cevfp if there is no use of the factor
    #printstyled("Removing missing factor payments\n", color=:yellow)
    delete.(mc.model, Array(qfe)[.!δ_evfp])
    delete.(mc.model, Array(qes)[.!δ_evfp])
    delete.(mc.model, Array(pes)[.!δ_evfp])
    delete.(mc.model, Array(pfe)[.!δ_evfp])
    delete.(mc.model, Array(peb)[.!δ_evfp])
    delete.(mc.model, Array(α_qes2[endws, :, :])[.!δ_evfp[endws, :, :]])
    delete.(mc.model, Array(α_qfe)[.!δ_evfp])
    delete.(mc.model, Array(qesf[endwf, acts, reg])[.!δ_evfp[endwf, acts, reg]])
    #delete.(mc.model, Array(evos)[.!δ_evfp])
    #delete.(mc.model, Array(evfp)[.!δ_evfp])
    if calibration
        delete.(mc.model, Array(σ_vff)[.!δ_evfp])
    end

    # Remove maks when there is no output produced
    #printstyled("Removing missing cross output\n", color=:yellow)
    #delete.(mc.model, Array(maks)[.!δ_maks])
    delete.(mc.model, Array(ps)[.!δ_maks])
    delete.(mc.model, Array(qca)[.!δ_maks])
    delete.(mc.model, Array(pca)[.!δ_maks])
    delete.(mc.model, Array(α_qca)[.!δ_maks])

    # Remove trade when no trade is allowed
    #printstyled("Removing missing trade\n", color=:yellow)
    delete.(mc.model, Array(α_qxs)[.!δ_qxs])
    delete.(mc.model, Array(qxs)[.!δ_qxs])
    delete.(mc.model, Array(pmds)[.!δ_qxs])
    delete.(mc.model, Array(pcif)[.!δ_qxs])
    delete.(mc.model, Array(pfob)[.!δ_qxs])
    delete.(mc.model, Array(ptrans)[.!δ_qxs.||.!δ_vtwr_sum])
    delete.(mc.model, Array(tms)[.!δ_qxs])
    delete.(mc.model, Array(txs)[.!δ_qxs])
    #delete.(mc.model, Array(vcif)[.!δ_qxs])
    #delete.(mc.model, Array(vfob)[.!δ_qxs])
    if calibration
        delete.(mc.model, Array(σ_qxs)[.!δ_qxs])
    end

    # Remove margins when no margins are allowed
    #printstyled("Removing missing margins\n", color=:yellow)
    delete.(mc.model, Array(qtmfsd)[.!δ_vtwr])
    if calibration
        delete.(mc.model, Array(σ_vtwr)[.!δ_vtwr])
    end
    delete.(mc.model, Array(α_qtmfsd)[.!δ_vtwr])
    #delete.(mc.model, Array(vtwr)[.!δ_vtwr])

    # Remove government consumption if not present
    #printstyled("Removing missing government spending\n", color=:yellow)
    delete.(mc.model, Array(qga)[.!δ_qga])
    delete.(mc.model, Array(pga)[.!δ_qga])
    delete.(mc.model, Array(α_qga)[.!δ_qga])
    delete.(mc.model, Array(qgm)[.!δ_qga])
    delete.(mc.model, Array(pgm)[.!δ_qga])
    delete.(mc.model, Array(qgd)[.!δ_qga])
    delete.(mc.model, Array(pgd)[.!δ_qga])
    delete.(mc.model, Array(α_qgdqgm["dom", :, :])[.!δ_qga])
    delete.(mc.model, Array(α_qgdqgm["imp", :, :])[.!δ_qga])
    delete.(mc.model, Array(γ_qgdqgm[:, :])[.!δ_qga])
    if calibration
        delete.(mc.model, Array(σ_vdg)[.!δ_qga])
        delete.(mc.model, Array(σ_vg)[.!δ_qga])
        delete.(mc.model, Array(ϵ_qgdqgm)[.!δ_qga])
    end
    delete.(mc.model, Array(tgm)[.!δ_qga])
    delete.(mc.model, Array(tgd)[.!δ_qga])
    #delete.(mc.model, Array(vdgp)[.!δ_qga])
    #delete.(mc.model, Array(vmgp)[.!δ_qga])

    # Remove investment consumption if not present
    #printstyled("Removing missing investment spending\n", color=:yellow)
    delete.(mc.model, Array(qia)[.!δ_qia])
    delete.(mc.model, Array(pia)[.!δ_qia])
    delete.(mc.model, Array(α_qia)[.!δ_qia])
    delete.(mc.model, Array(qim)[.!δ_qia])
    delete.(mc.model, Array(pim)[.!δ_qia])
    delete.(mc.model, Array(qid)[.!δ_qia])
    delete.(mc.model, Array(pid)[.!δ_qia])
    delete.(mc.model, Array(α_qidqim["dom", :, :])[.!δ_qia])
    delete.(mc.model, Array(α_qidqim["imp", :, :])[.!δ_qia])
    delete.(mc.model, Array(γ_qidqim[:, :])[.!δ_qia])
    if calibration
        delete.(mc.model, Array(σ_vdi)[.!δ_qia])
        delete.(mc.model, Array(σ_vi)[.!δ_qia])
        delete.(mc.model, Array(ϵ_qidqim)[.!δ_qia])
    end
    delete.(mc.model, Array(tim)[.!δ_qia])
    delete.(mc.model, Array(tid)[.!δ_qia])
    #delete.(mc.model, Array(vdip)[.!δ_qia])
    #delete.(mc.model, Array(vmip)[.!δ_qia])


    # Remove investment consumption if not present
    #printstyled("Removing missing firm spending\n", color=:yellow)
    delete.(mc.model, Array(qfa)[.!δ_qfa])
    delete.(mc.model, Array(pfa)[.!δ_qfa])
    delete.(mc.model, Array(α_qfa)[.!δ_qfa])
    delete.(mc.model, Array(qfm)[.!δ_qfa])
    delete.(mc.model, Array(pfm)[.!δ_qfa])
    delete.(mc.model, Array(qfd)[.!δ_qfa])
    delete.(mc.model, Array(pfd)[.!δ_qfa])
    delete.(mc.model, Array(α_qfdqfm["dom", :, :, :])[.!δ_qfa])
    delete.(mc.model, Array(α_qfdqfm["imp", :, :, :])[.!δ_qfa])
    delete.(mc.model, Array(γ_qfdqfm[:, :, :])[.!δ_qfa])
    if calibration
        delete.(mc.model, Array(σ_vdf)[.!δ_qfa])
        delete.(mc.model, Array(σ_vf)[.!δ_qfa])
        delete.(mc.model, Array(ϵ_qfdqfm)[.!δ_qfa])
    end
    delete.(mc.model, Array(tfm)[.!δ_qfa])
    delete.(mc.model, Array(tfd)[.!δ_qfa])
    #delete.(mc.model, Array(vdfp)[.!δ_qfa])
    #delete.(mc.model, Array(vmfp)[.!δ_qfa])
    # All model equations
    @constraints(mc.model,
        begin
            # Firms (top nest)
            e_qintva[a=acts, r=reg; δ_act[a, r]], log.([qint[a, r], qva[a, r]]) .== log.(ces(qo[a, r], [pint[a, r], pva[a, r]], α_qintva[:, a, r], esubt[a, r], γ_qintva[a, r]))
            e_qo[a=acts, r=reg; δ_act[a, r]], log.(qo[a, r] * po[a, r]) == log(qva[a, r] * pva[a, r] + qint[a, r] * pint[a, r])

            # Firms (second nest)
            e_qfa[a=acts, r=reg; δ_act[a, r]], log.(Vector(qfa[:, a, r])[δ_qfa[:, a, r]]) .== log.(ces(qint[a, r], Vector(pfa[:, a, r])[δ_qfa[:, a, r]], Vector(α_qfa[:, a, r])[δ_qfa[:, a, r]], esubc[a, r], γ_qfa[a, r]))
            e_pint[a=acts, r=reg; δ_act[a, r]], log.(qint[a, r] * pint[a, r]) == log.(sum(Vector(pfa[:, a, r] .* qfa[:, a, r])[δ_qfa[:, a, r]]))
            e_qfe[a=acts, r=reg; δ_act[a, r]], log.(Vector(qfe[:, a, r])[δ_evfp[:, a, r]]) .== log.(ces(qva[a, r], Vector(pfe[:, a, r])[δ_evfp[:, a, r]], Vector(α_qfe[:, a, r])[δ_evfp[:, a, r]], esubva[a, r], γ_qfe[a, r]))

            e_pva[a=acts, r=reg; δ_act[a, r]], log.(qva[a, r] * pva[a, r]) == log.(sum(Vector(pfe[:, a, r] .* qfe[:, a, r])[δ_evfp[:, a, r]]))
            e_qfdqfm[c=comm, a=acts, r=reg; δ_qfa[c, a, r] && δ_act[a, r]], log.([qfd[c, a, r], qfm[c, a, r]]) .== log.(ces(qfa[c, a, r], [pfd[c, a, r], pfm[c, a, r]], α_qfdqfm[:, c, a, r], esubd[c, r], γ_qfdqfm[c, a, r]))
            e_pfa[c=comm, a=acts, r=reg; δ_qfa[c, a, r] && δ_act[a, r]], log(pfa[c, a, r] * qfa[c, a, r]) == log(qfd[c, a, r] * pfd[c, a, r] + qfm[c, a, r] * pfm[c, a, r])

            # Firms (distribution)
            e_qca[a=acts, r=reg; δ_act[a, r]], log.(Vector(qca[:, a, r])[δ_maks[:, a, r]]) .== log.(Vector(ces(qo[a, r], Vector(ps[:, a, r])[δ_maks[:, a, r]], Vector(α_qca[:, a, r])[δ_maks[:, a, r]], etraq[a, r], γ_qca[a, r])))
            e_po[a=acts, r=reg; δ_act[a, r]], log.(po[a, r] * qo[a, r]) == log.(sum(Vector(qca[:, a, r] .* ps[:, a, r])[δ_maks[:, a, r]]))
            e_pca[c=comm, r=reg], log.((esubq[c, r] == 0 ? Vector(pca[c, :, r])[δ_maks[c, :, r]] : Vector(qca[c, :, r])[δ_maks[c, :, r]])) .== log.((esubq[c, r] == 0 ? pds[c, r] : Vector(ces(qc[c, r], Vector(pca[c, :, r])[δ_maks[c, :, r]], Vector(α_pca[c, :, r])[δ_maks[c, :, r]], 1 / esubq[c, r], γ_pca[c, r]))))
            e_qc[c=comm, r=reg], log(pds[c, r] * qc[c, r]) == log(sum(Vector(pca[c, :, r] .* qca[c, :, r])[δ_maks[c, :, r]]))
            e_ps[c=comm, a=acts, r=reg; δ_maks[c, a, r]], log(pca[c, a, r]) == log(ps[c, a, r] * to[c, a, r])

            # Endowments
            e_peb[e=endw, a=acts, r=reg; δ_evfp[e, a, r]], log(qfe[e, a, r]) == log(qes[e, a, r])
            e_pfe[e=endw, a=acts, r=reg; δ_evfp[e, a, r]], log(pfe[e, a, r]) == log(peb[e, a, r] * tfe[e, a, r])
            e_pfactor[r=reg], log(pfactor[r] * sum(Array(qfe[:, :, r])[δ_evfp[:, :, r]])) == log(sum(Array(qfe[:, :, r] .* peb[:, :, r])[δ_evfp[:, :, r]]))

            # Income
            e_fincome[r=reg], log(fincome[r]) == log(sum(Array(peb[:, :, r] .* qes[:, :, r])[δ_evfp[:, :, r]]) - δ[r] * pinv[r] * kb[r])
            e_y[r=reg], log(y[r]) ==
                        log(
                fincome[r] +
                sum(qpd[:, r] .* pds[:, r] .* (tpd[:, r] .- 1)) +
                sum(qpm[:, r] .* pms[:, r] .* (tpm[:, r] .- 1)) +
                sum(Vector(qgd[:, r] .* pds[:, r] .* (tgd[:, r] .- 1))[δ_qga[:, r]]) +
                sum(Vector(qgm[:, r] .* pms[:, r] .* (tgm[:, r] .- 1))[δ_qga[:, r]]) +
                sum(Vector(qid[:, r] .* pds[:, r] .* (tid[:, r] .- 1))[δ_qia[:, r]]) +
                sum(Vector(qim[:, r] .* pms[:, r] .* (tim[:, r] .- 1))[δ_qia[:, r]]) +
                sum(Array(qfd[:, :, r] .* pfd[:, :, r] ./ tfd[:, :, r] .* (tfd[:, :, r] .- 1))[δ_qfa[:, :, r]]) +
                sum(Array(qfm[:, :, r] .* pfm[:, :, r] ./ tfm[:, :, r] .* (tfm[:, :, r] .- 1))[δ_qfa[:, :, r]]) +
                sum(Array(qca[:, :, r] .* ps[:, :, r] .* (to[:, :, r] .- 1))[δ_maks[:, :, r]]) +
                sum(Array(qfe[:, :, r] .* peb[:, :, r] .* (tfe[:, :, r] .- 1))[δ_evfp[:, :, r]]) +
                sum(Array(qxs[:, r, :] .* pfob[:, r, :] ./ txs[:, r, :] .* (txs[:, r, :] .- 1))[δ_qxs[:, r, :]]) +
                sum(Array(qxs[:, :, r] .* pcif[:, :, r] .* (tms[:, :, r] .- 1))[δ_qxs[:, :, r]])
            )

            # Utility
            e_ug, log.(ug) .== log.(yg ./ pop ./ pgov)
            e_p, log.(p) .== log.(ppriv .* σyp .+ pgov .* σyg .+ psave .* (1 .- σyp .- σyg))
            e_u, log.(u) .== log.((y ./ p ./ pop)) 

            # Household Income
            e_yp, log.(yp) .== log.(y .* Vector(σyp) .* uelas ./ uepriv)

            # Household consumption
            e_qpa[r=reg], log.([Vector(qpa[:, r] ./ pop[r]); 1]) .== log.(cde(Vector(1 .- subpar[:, r]), Vector(β_qpa[:, r]), Vector(incpar[:, r]), up[r], Vector(ppa[:, r]), yp[r] / pop[r]))
            e_uepriv[r=reg], uepriv[r] == sum(qpa[:, r] .* ppa[:, r] .* Vector(incpar[:, r])) / yp[r]
            e_uelas[r=reg], uelas[r] == 1 / (σyp[r] / uepriv[r] + σyg[r] + (1 - σyp[r] - σyg[r]))
            e_qpdqpm[c=comm, r=reg], log.([qpd[c, r], qpm[c, r]]) .== log.(ces(qpa[c, r], [ppd[c, r], ppm[c, r]], α_qpdqpm[:, c, r], esubd[c, r], γ_qpdqpm[c, r]))
            e_ppa, log.(qpa .* ppa) .== log.(ppd .* qpd .+ ppm .* qpm)
            e_ppriv[r=reg], log(ppriv[r] * sum(qpa[:, r])) == log(sum(ppa[:, r] .* qpa[:, r]))

            # Government Income
            e_yg, log.(yg) .== log.(y .* σyg .* uelas)

            # Government expenditure
            e_qga[c=comm, r=reg; δ_qga[c, r]], log(pga[c, r] * qga[c, r]) == log(yg[r] * α_qga[c, r]) ##This one
            e_pgov[r=reg], log(pgov[r] * sum(Vector(qga[:, r])[δ_qga[:, r]])) == log(sum(Vector(qga[:, r] .* pga[:, r])[δ_qga[:, r]]))
            e_qgdqgm[c=comm, r=reg; δ_qga[c, r]], log.([qgd[c, r], qgm[c, r]]) .== log.(ces(qga[c, r], [pgd[c, r], pgm[c, r]], α_qgdqgm[:, c, r], esubd[c, r], γ_qgdqgm[c, r]))
            e_pga[c=comm, r=reg; δ_qga[c, r]], log.(qga[c, r] .* pga[c, r]) .== log.(pgd[c, r] .* qgd[c, r] .+ pgm[c, r] .* qgm[c, r])

            # Saving
            e_qsave, log.(y .* uelas) .== log.(σyp .* y .* uelas .+ σyg .* y .* uelas .+ psave .* qsave)

            # Investment consumption
            e_qia[r=reg], log.(Vector(qia[:, r])[δ_qia[:, r]]) .== log.(ces(qinv[r], Vector(pia[:, r])[δ_qia[:, r]], Vector(α_qia[:, r])[δ_qia[:, r]], 0, γ_qia[r]))
            e_pinv[r=reg], log.(pinv[r] * sum(Vector(qia[:, r])[δ_qia[:, r]])) == log.(sum(Vector(pia[:, r] .* qia[:, r])[δ_qia[:, r]]))
            e_qidqim[c=comm, r=reg; δ_qia[c, r]], log.([qid[c, r], qim[c, r]]) .== log.(ces(qia[c, r], [pid[c, r], pim[c, r]], Vector(α_qidqim[:, c, r]), esubd[c, r], γ_qidqim[c, r]))
            e_pia[c=comm, r=reg; δ_qia[c, r]], log(pia[c, r] .* qia[c, r]) == log(pid[c, r] * qid[c, r] + pim[c, r] * qim[c, r])
            #e_psave[r=reg], log(psave[r]) == log(pinv[r] + sum(((qinv .- δ .* kb) .- qsave) .* pinv ./ sum(qinv .- δ .* kb)))
            e_psave[r=reg], log(psave[r]) == log(pinv[r]) + sum(((qinv .- δ .* kb) .- qsave) .* log.(pinv) ./ sum(qinv .- δ .* kb))

            # Trade - exports
            e_qms[c=comm, r=reg], log.(qms[c, r]) == log.(sum(Vector(qfm[c, :, r])[δ_qfa[c, :, r]]) + qpm[c, r] + (δ_qga[c, r] ? qgm[c, r] : 0) + (δ_qia[c, r] ? qim[c, r] : 0))
            e_qxs[c=comm, r=reg], log.(Array(qxs[c, :, r])[δ_qxs[c, :, r]]) .== log.(ces(qms[c, r], Vector(pmds[c, :, r])[δ_qxs[c, :, r]], Vector(α_qxs[c, :, r])[δ_qxs[c, :, r]], esubm[c, r], γ_qxs[c, r]))
            e_pms[c=comm, r=reg], log(pms[c, r] * qms[c, r]) == log(sum(Vector(pmds[c, :, r] .* qxs[c, :, r])[δ_qxs[c, :, r]]))

            # Trade - margins
            e_qtmfsd[m=marg, c=comm, s=reg, d=reg; δ_vtwr[m, c, s, d]], log(qtmfsd[m, c, s, d]) == log(α_qtmfsd[m, c, s, d] * qxs[c, s, d])
            e_ptrans[c=comm, s=reg, d=reg; δ_vtwr_sum[c, s, d]], log(ptrans[c, s, d] * sum(Array(qtmfsd[:, c, s, d])[δ_vtwr[:, c, s, d]])) == log(sum(Array(qtmfsd[:, c, s, d] .* pt[:])[δ_vtwr[:, c, s, d]]))
            e_qtm[m=marg], log(qtm[m]) == log(sum(Array(qtmfsd[m, :, :, :])[δ_vtwr[m, :, :, :]]))
            e_qst[m=marg], log.(qst[m, :]) .== log.(ces(qtm[m], pds[m, :], Vector(α_qst[m, :]), esubs[m], γ_qst[m]))
            e_pt[m=marg], log.(pt[m] * qtm[m]) == log.(sum(pds[m, :] .* qst[m, :]))

            # Trade - imports / exports
            e_pfob[c=comm, s=reg, d=reg; δ_qxs[c, s, d]], log(pfob[c, s, d]) == log(pds[c, s] * tx[c, s] * txs[c, s, d])
            e_pcif[c=comm, s=reg, d=reg; δ_qxs[c, s, d]], log(pcif[c, s, d] * qxs[c, s, d]) == log(pfob[c, s, d] * (qxs[c, s, d]) + (δ_vtwr_sum[c, s, d] ? ptrans[c, s, d] * sum(Vector(qtmfsd[:, c, s, d])[δ_vtwr[:, c, s, d]]) : 0))
            e_pmds[c=comm, s=reg, d=reg; δ_qxs[c, s, d]], log(pmds[c, s, d]) == log(pcif[c, s, d] * tm[c, d] * tms[c, s, d])

            # Domestic market clearing 
            e_qds[c=comm, r=reg], log(qds[c, r]) == log(sum(Vector(qfd[c, :, r])[δ_qfa[c, :, r]]) + qpd[c, r] + (δ_qga[c, r] ? qgd[c, r] : 0) + (δ_qia[c, r] ? qid[c, r] : 0))
            e_pds[c=comm, r=reg], log(qc[c, r]) == log(qds[c, r] + sum(Vector(qxs[c, r, :])[δ_qxs[c, r, :]]) + (c ∈ marg ? qst[c, r] : 0))

            # Taxes
            e_pfd[c=comm, a=acts, r=reg; δ_qfa[c, a, r]], log(pfd[c, a, r]) == log(pds[c, r] * tfd[c, a, r])
            e_pfm[c=comm, a=acts, r=reg; δ_qfa[c, a, r]], log(pfm[c, a, r]) == log(pms[c, r] * tfm[c, a, r])
            e_ppd, log.(ppd) .== log.(pds .* tpd)
            e_ppm, log.(ppm) .== log.(pms .* tpm)
            e_pgd[c=comm, r=reg; δ_qga[c, r]], log(pgd[c, r]) == log(pds[c, r] * tgd[c, r])
            e_pgm[c=comm, r=reg; δ_qga[c, r]], log(pgm[c, r]) == log(pms[c, r] * tgm[c, r])
            e_pid[c=comm, r=reg; δ_qia[c, r]], log(pid[c, r]) == log(pds[c, r] * tid[c, r])
            e_pim[c=comm, r=reg; δ_qia[c, r]], log(pim[c, r]) == log(pms[c, r] * tim[c, r])

            # Factor Market
            e_pe1[e=endwm, r=reg], log.(qe[e, r]) == log.(sum(Vector(qfe[e, :, r])[δ_evfp[e, :, r]]))
            e_qes1[e=endwm, a=acts, r=reg; δ_evfp[e, a, r]], log(pes[e, a, r]) == log(pe[e, r])
            e_qes2[e=endws, r=reg], log.(Vector(qes[e, :, r])[δ_evfp[e, :, r]]) .== log.(Vector(ces(qe[e, r], Vector(pes[e, :, r])[δ_evfp[e, :, r]], Vector(α_qes2[e, :, r])[δ_evfp[e, :, r]], etrae[e, r], γ_qes2[e, r])))
            e_pe2[e=endws, r=reg], log(pe[e, r] * qe[e, r]) == log(sum(Vector(pes[e, :, r] .* qes[e, :, r])[δ_evfp[e, :, r]]))
            e_qes3[e=endwf, a=acts, r=reg; δ_evfp[e, a, r]], log(qes[e, a, r]) == log(qesf[e, a, r])
            e_pes[e=endw, a=acts, r=reg; δ_evfp[e, a, r]], log(peb[e, a, r]) == log(pes[e, a, r] * tinc[e, a, r])

            # Investment is a fixed share of global investment
            e_pcgdswld, log(pcgdswld) == log(sum(pinv .* (qinv .- δ .* kb)) / sum(qinv .- δ .* kb))
            e_walras_sup, log(walras_sup) == log(pcgdswld * globalcgds)
            e_walras_dem, log(walras_dem) == log(sum(psave .* qsave))

            # Pfactwld
            e_rorg, log(pfactwld * sum(Array(qfe)[δ_evfp])) == log(sum(Array(peb .* qfe)[δ_evfp]))

            # Capital accumulation
            e_kb[r=reg], log(ρ[r] * kb[r]) == log(sum(qe[endwc, r]))
            e_ke, log.(ke) .== log.(qinv .+ (1 .- δ) .* kb)

            e_rore[r=reg], log(α_qinv[r] * rore[r]) == log(rorc[r] / (ke[r] / kb[r])^rorflex[r])
            e_rorc[r=reg], log(rorc[r] * (sum(Array(qes[endwc, :, r] )[δ_evfp[endwc, :, r]]) - δ[r] .* kb[r])) == log(sum(Array(qes[endwc, :, r] )[δ_evfp[endwc, :, r]]) * (rental[r] / pinv[r]))
            e_rental[r=reg], log(sum(Array(qes[endwc, :, r] )[δ_evfp[endwc, :, r]])*rental[r]) == log(sum(Array(qes[endwc, :, r] .* pes[endwc, :, r])[δ_evfp[endwc, :, r]]))
        end
    )

    if rordelta == 1
        @constraints(mc.model,
            begin
                #e_rorc[r=reg], log((sum(Array(pes[endwc, :, r] .* qes[endwc, :, r])[δ_evfp[endwc, :, r]]) - δ[r] * pinv[r] * kb[r]) * rorc[r]) == log(sum(Array(pes[endwc, :, r] .* qes[endwc, :, r])[δ_evfp[endwc, :, r]]) * rental[r] / pinv[r])
                e_qinv[r=reg], log(rore[r]) == log(rorg)
                e_globalcgds, log(globalcgds) == log(sum(qinv .- δ .* kb))
            end
        )
    else
        @constraints(mc.model,
            begin
                e_qinv, log.(qinv) .== log.(Vector(α_qinv) .* globalcgds .+ δ .* kb)
                e_globalcgds, log(rorg * sum(qinv .- δ .* kb)) == log(sum((qinv .- δ .* kb) .* rore))
            end
        )
    end

    if calibration

        @constraints(mc.model,
            begin

                # Values
                cvdpp, log.(vdpp) .== log.(ppd .* qpd)

                # Soft parameter constraints
                sf_α_qxs[c=comm, d=reg], log(sum(Vector(α_qxs[c, :, d])[δ_qxs[c, :, d]])) == log(ϵ_qxs[c, d])
                sf_α_qfe[a=acts, r=reg; δ_act[a, r]], log(sum(Vector(α_qfe[:, a, r])[δ_evfp[:, a, r]])) == log(ϵ_qfe[a, r])
                sf_α_qes2[e=endws, r=reg], log(sum(Vector(α_qes2[e, :, r])[δ_evfp[e, :, r]])) == log(ϵ_qes2[e, r])
                sf_α_qfdqfm[c=comm, a=acts, r=reg; δ_qfa[c, a, r]], log(sum(α_qfdqfm[:, c, a, r])) == log(ϵ_qfdqfm[c, a, r])
                sf_α_qpdqpm[c=comm, r=reg], log(sum(α_qpdqpm[:, c, r])) == log(ϵ_qpdqpm[c, r])
                sf_α_qgdqgm[c=comm, r=reg; δ_qga[c, r]], log(sum(α_qgdqgm[:, c, r])) == log(ϵ_qgdqgm[c, r])
                sf_α_qidqim[c=comm, r=reg; δ_qia[c, r]], log(sum(α_qidqim[:, c, r])) == log(ϵ_qidqim[c, r])
                sf_α_qia[r=reg], log(sum(Vector(α_qia[:, r])[δ_qia[:, r]])) == log(ϵ_qia[r])
                sf_α_qga[r=reg], log(sum(Vector(α_qga[:, r])[δ_qga[:, r]])) == log(ϵ_qga[r])
                sf_α_qfa[a=acts, r=reg; δ_act[a,r]], log(sum(Vector(α_qfa[:, a, r])[δ_qfa[:, a, r]])) == log(ϵ_qfa[a, r])
                sf_α_qintva[a=acts, r=reg; δ_act[a,r]], log(sum(α_qintva[["int", "va"], a, r])) == log(ϵ_qintva[a, r])
                sf_α_qinv, log(sum(α_qinv)) == log(ϵ_qinv)
                sf_save, log.(σsave .+ σyp .+ σyg) .== log(1)

                # Shares (helpers)
                e_σ_vp[c=comm, r=reg], log(σ_vp[c, r]) + log(sum(ppd[:, r] .* qpd[:, r] .+ ppm[:, r] .* qpm[:, r])) == log(ppd[c, r] .* qpd[c, r] + ppm[c, r] .* qpm[c, r])
                e_σ_vdp[c=comm, r=reg], log(σ_vdp[c, r]) + log(ppd[c, r] .* qpd[c, r] .+ ppm[c, r] .* qpm[c, r]) == log(ppd[c, r] .* qpd[c, r])
                e_σ_vg[c=comm, r=reg; δ_qga[c, r]], log(σ_vg[c, r]) + log(sum(Vector(pgd[:, r] .* qgd[:, r] .+ pgm[:, r] .* qgm[:, r])[δ_qga[:, r]])) == log(pgd[c, r] * qgd[c, r] + pgm[c, r] * qgm[c, r])
                e_σ_vdg[c=comm, r=reg; δ_qga[c, r]], log(σ_vdg[c, r]) + log(pgd[c, r] .* qgd[c, r] .+ pgm[c, r] .* qgm[c, r]) == log(pgd[c, r] .* qgd[c, r])
                e_σ_vi[c=comm, r=reg; δ_qia[c, r]], log(σ_vi[c, r]) + log(sum(Vector(pid[:, r] .* qid[:, r] .+ pim[:, r] .* qim[:, r])[δ_qia[:, r]])) == log(pid[c, r] * qid[c, r] + pim[c, r] * qim[c, r])
                e_σ_vdi[c=comm, r=reg; δ_qia[c, r]], log(σ_vdi[c, r]) + log(pid[c, r] .* qid[c, r] .+ pim[c, r] .* qim[c, r]) == log(pid[c, r] .* qid[c, r])
                e_σ_vf[c=comm, a=acts, r=reg; δ_qfa[c, a, r]], log(σ_vf[c, a, r]) + log(sum(Vector(pfd[:, a, r] .* qfd[:, a, r] .+ pfm[:, a, r] .* qfm[:, a, r])[δ_qfa[:, a, r]])) == log(pfd[c, a, r] .* qfd[c, a, r] + pfm[c, a, r] .* qfm[c, a, r])
                e_σ_vdf[c=comm, a=acts, r=reg; δ_qfa[c, a, r]], log(σ_vdf[c, a, r]) + log(pfd[c, a, r] .* qfd[c, a, r] .+ pfm[c, a, r] .* qfm[c, a, r]) == log(pfd[c, a, r] .* qfd[c, a, r])
                e_σ_vif[a=acts, r=reg; δ_act[a,r]], log(σ_vif[a, r]) + log(pva[a, r] * qva[a, r] + pint[a, r] * qint[a, r]) == log(pint[a, r] * qint[a, r])
                e_σ_vff[e=endw, a=acts, r=reg; δ_evfp[e, a, r]], log(σ_vff[e, a, r]) + log(sum(Vector(pfe[:, a, r] .* qfe[:, a, r])[δ_evfp[:, a, r]])) == log(pfe[e, a, r] * qfe[e, a, r])
                e_σ_vtwr[m=marg, c=comm, s=reg, d=reg; δ_vtwr[m, c, s, d]], log(σ_vtwr[m, c, s, d]) + log(pcif[c, s, d] * qxs[c, s, d]) == log(pt[m] * qtmfsd[m, c, s, d])
                e_σ_qxs[c=comm, s=reg, d=reg; δ_qxs[c, s, d]], log(σ_qxs[c, s, d]) + log(sum(Vector(pcif[c, :, d] .* qxs[c, :, d])[δ_qxs[c, :, d]])) == log(pcif[c, s, d] .* qxs[c, s, d])

                e_σ_qinv[r=reg], σ_qinv[r] == psave[r] * qsave[r] / sum(psave .* qsave)
                e_σ_ρ[r=reg], log((kb[r] * pinv[r]) * σ_ρ[r]) == log(sum(Array(pes[:, :, r] .* qes[:, :, r])[δ_evfp[:, :, r]]))
            end
        )

    end

    set_attribute(mc.model, "max_iter", max_iter)
    set_attribute(mc.model, "constr_viol_tol", constr_viol_tol)
    set_attribute(mc.model, "bound_push", bound_push)

    return nothing
end
