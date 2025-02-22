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
function build_model!(mc; max_iter=50, constr_viol_tol=1e-8, bound_push=1e-15)
    # Structural parameters (some CES/CET options are not happening)
    δ_evfp = .!isnan.(mc.data["α_qfe"]) .&& mc.data["α_qfe"] .!= 0
    δ_maks = .!isnan.(mc.data["α_qca"]) .&& mc.data["α_qca"] .!= 0
    δ_vtwr = .!isnan.(mc.data["α_qtmfsd"]) .&& mc.data["α_qtmfsd"] .!= 0
    δ_vtwr_sum = NamedArray(.!mapslices(all, isnan.(mc.data["α_qtmfsd"]) .|| mc.data["α_qtmfsd"] .== 0, dims=1)[1, :, :, :], names(mc.data["vtwr"])[[2, 3, 4]])
    δ_qxs = .!isnan.(mc.data["α_qxs"]) .&& mc.data["α_qxs"] .!= 0
    δ_qga = .!isnan.(mc.data["α_qga"]) .&& mc.data["α_qga"] .!= 0
    δ_qia = .!isnan.(mc.data["α_qia"]) .&& mc.data["α_qia"] .!= 0

    # Read  sets
    (; reg, comm, marg, acts, endw, endwc, endws, endwm, endwms, endwf) = NamedTuple(Dict(Symbol(k) => mc.sets[k] for k ∈ keys(mc.sets)))

    # Read hard parameters
    (; esubt, esubc, esubva, esubd, etraq, esubq, subpar, incpar, etrae, esubg, esubm, esubs) = NamedTuple(Dict(Symbol(k) => mc.parameters[k] for k ∈ keys(mc.parameters)))

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
            u[reg]
            ug[reg]
            us[reg]

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

            # Soft parameters
            α_qintva[["int", "va"], acts, reg]
            ϵ_qintva[acts, reg]
            γ_qintva[acts, reg]
            α_qfa[comm, acts, reg]
            ϵ_qfa[acts, reg]
            γ_qfa[acts, reg]
            α_qfe[endw, acts, reg]
            γ_qfe[acts, reg]
            ϵ_qfe[acts, reg]
            α_qfdqfm[["dom", "imp"], comm, acts, reg]
            γ_qfdqfm[comm, acts, reg]
            ϵ_qfdqfm[comm, acts, reg]
            α_qca[comm, acts, reg]
            γ_qca[acts, reg]
            α_pca[comm, acts, reg]
            γ_pca[acts, reg]
            σyp[reg]
            σyg[reg]
            σsave[reg]
            α_qga[comm, reg]
            ϵ_qga[reg]
            γ_qga[reg]
            α_qia[comm, reg]
            γ_qia[reg]
            ϵ_qia[reg]
            α_qpdqpm[["dom", "imp"], acts, reg]
            γ_qpdqpm[acts, reg]
            ϵ_qpdqpm[acts, reg]
            α_qgdqgm[["dom", "imp"], acts, reg]
            γ_qgdqgm[acts, reg]
            ϵ_qgdqgm[acts, reg]
            β_qpa[comm, reg]
            α_qidqim[["dom", "imp"], acts, reg]
            γ_qidqim[acts, reg]
            ϵ_qidqim[acts, reg]
            α_qxs[comm, reg, reg]
            γ_qxs[comm, reg]
            ϵ_qxs[comm, reg]
            α_qtmfsd[marg, comm, reg, reg]
            α_qst[marg, reg]
            γ_qst[marg]
            α_qes2[endws, acts, reg]
            γ_qes2[endws, reg]
            ϵ_qes2[endws, reg]
            α_qinv[reg]
            ϵ_qinv

            δ[reg]
            ρ[reg]

            # Values
            vdfp[comm, acts, reg]
            vmfp[comm, acts, reg]
            vdpp[comm, reg]
            vmpp[comm, reg]
            vdgp[comm, reg]
            vmgp[comm, reg]
            vdip[comm, reg]
            vmip[comm, reg]
            evfp[endw, acts, reg]
            evos[endw, acts, reg]
            vfob[comm, reg, reg]
            vcif[comm, reg, reg]
            vst[marg, reg]
            vtwr[marg, comm, reg, reg]
            maks[comm, acts, reg]
            vkb[reg]

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


    # All model equations
    @constraints(mc.model,
        begin
            # Firms (top nest)
            e_qintva[a=acts, r=reg], log.([qint[a, r], qva[a, r]]) .== log.(ces(qo[a, r], [pint[a, r], pva[a, r]], α_qintva[:, a, r], esubt[a, r], γ_qintva[a, r]))
            e_qo, log.(qo .* po) .== log.(qva .* pva .+ qint .* pint)

            # Firms (second nest)
            e_qfa[a=acts, r=reg], log.(qfa[:, a, r]) .== log.(ces(qint[a, r], pfa[:, a, r], Vector(α_qfa[:, a, r]), esubc[a, r], γ_qfa[a, r]))
            e_pint[a=acts, r=reg], log.(qint[a, r] * pint[a, r]) == log.(sum(pfa[:, a, r] .* qfa[:, a, r]))
            e_qfe[a=acts, r=reg], log.(Vector(qfe[:, a, r])[δ_evfp[:, a, r]]) .== log.(ces(qva[a, r], Vector(pfe[:, a, r])[δ_evfp[:, a, r]], Vector(α_qfe[:, a, r])[δ_evfp[:, a, r]], esubva[a, r], γ_qfe[a, r]))

            e_pva[a=acts, r=reg], log.(qva[a, r] * pva[a, r]) == log.(sum(pfe[:, a, r] .* qfe[:, a, r]))
            e_qfdqfm[c=comm, a=acts, r=reg], log.([qfd[c, a, r], qfm[c, a, r]]) .== log.(ces(qfa[c, a, r], [pfd[c, a, r], pfm[c, a, r]], α_qfdqfm[:, c, a, r], esubd[c, r], γ_qfdqfm[c, a, r]))
            e_pfa, log.(pfa .* qfa) .== log.(qfd .* pfd .+ qfm .* pfm)

            # Firms (distribution)
            e_qca[a=acts, r=reg], log.(Vector(qca[:, a, r])[δ_maks[:, a, r]]) .== log.(Vector(ces(qo[a, r], Vector(ps[:, a, r])[δ_maks[:, a, r]], Vector(α_qca[:, a, r])[δ_maks[:, a, r]], etraq[a, r], γ_qca[a, r])))
            e_po[a=acts, r=reg], log.(po[a, r] * qo[a, r]) == log.(sum(qca[:, a, r] .* ps[:, a, r]))
            e_pca[c=comm, r=reg], log.((esubq[c, r] == 0 ? Vector(pca[c, :, r])[δ_maks[c, :, r]] : Vector(qca[c, :, r])[δ_maks[c, :, r]])) .== log.((esubq[c, r] == 0 ? pds[c, r] : Vector(ces(qc[c, r], Vector(pca[c, :, r])[δ_maks[c, :, r]], Vector(α_pca[c, :, r])[δ_maks[c, :, r]], 1 / esubq[c, r], γ_pca[c, r]))))
            e_qc[c=comm, r=reg], log(pds[c, r] * qc[c, r]) == log(sum(pca[c, :, r] .* qca[c, :, r]))
            e_ps[c=comm, a=acts, r=reg; δ_maks[c, a, r]], log(pca[c, a, r]) == log(ps[c, a, r] * to[c, a, r])

            # Endowments
            e_peb[e=endw, a=acts, r=reg; δ_evfp[e, a, r]], log(qfe[e, a, r]) == log(qes[e, a, r])
            e_pfe[e=endw, a=acts, r=reg; δ_evfp[e, a, r]], log(pfe[e, a, r]) == log(peb[e, a, r] * tfe[e, a, r])
            e_pfactor[r=reg], log(pfactor[r] * sum(qfe[:, :, r])) == log(sum(qfe[:, :, r] .* peb[:, :, r]))

            # Income
            e_fincome[r=reg], log(fincome[r]) == log(sum(peb[:, :, r] .* qes[:, :, r]) .- δ[r] .* pinv[r] .* kb[r])
            e_y[r=reg], log(y[r]) ==
                        log(
                fincome[r] +
                sum(qpd[:, r] .* pds[:, r] .* (tpd[:, r] .- 1)) +
                sum(qpm[:, r] .* pms[:, r] .* (tpm[:, r] .- 1)) +
                sum(qgd[:, r] .* pds[:, r] .* (tgd[:, r] .- 1)) +
                sum(qgm[:, r] .* pms[:, r] .* (tgm[:, r] .- 1)) +
                sum(qid[:, r] .* pds[:, r] .* (tid[:, r] .- 1)) +
                sum(qim[:, r] .* pms[:, r] .* (tim[:, r] .- 1)) +
                sum(qfd[:, :, r] .* pfd[:, :, r] ./ tfd[:, :, r] .* (tfd[:, :, r] .- 1)) +
                sum(qfm[:, :, r] .* pfm[:, :, r] ./ tfm[:, :, r] .* (tfm[:, :, r] .- 1)) +
                sum(qca[:, :, r] .* ps[:, :, r] .* (to[:, :, r] .- 1)) +
                sum(qfe[:, :, r] .* peb[:, :, r] .* (tfe[:, :, r] .- 1)) +
                sum(qxs[:, r, :] .* pfob[:, r, :] ./ txs[:, r, :] .* (txs[:, r, :] .- 1)) +
                sum(qxs[:, :, r] .* pcif[:, :, r] .* (tms[:, :, r] .- 1))
            )


            # Utility
            e_ug, ug .== yg ./ pop ./ pgov
            e_us, us .== qsave ./ pop
            e_u, u .== up .^ σyp .* ug .^ σyg .* us .^ (1 .- σyp .- σyg)

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
            e_pgov[r=reg], log(pgov[r] * sum(qga[:, r])) == log(sum(qga[:, r] .* pga[:, r]))
            e_qgdqgm[c=comm, r=reg; δ_qga[c, r]], log.([qgd[c, r], qgm[c, r]]) .== log.(ces(qga[c, r], [pgd[c, r], pgm[c, r]], α_qgdqgm[:, c, r], esubd[c, r], γ_qgdqgm[c, r]))
            e_pga[c=comm, r=reg; δ_qga[c, r]], log.(qga[c, r] .* pga[c, r]) .== log.(pgd[c, r] .* qgd[c, r] .+ pgm[c, r] .* qgm[c, r])

            # Saving
            e_qsave, log.(y .* uelas) .== log.(σyp .* y .* uelas .+ σyg .* y .* uelas .+ psave .* qsave)

            # Investment consumption
            e_qia[r=reg], log.(Vector(qia[:, r])[δ_qia[:,r]]) .== log.(ces(qinv[r], Vector(pia[:, r])[δ_qia[:,r]], Vector(α_qia[:, r])[δ_qia[:,r]], 0, γ_qia[r]))
            e_pinv[r=reg], log.(pinv[r] * sum(qia[:, r])) == log.(sum(pia[:, r] .* qia[:, r]))
            e_qidqim[c=comm, r=reg;δ_qia[c,r]], log.([qid[c, r], qim[c, r]]) .== log.(ces(qia[c, r], [pid[c, r], pim[c, r]], Vector(α_qidqim[:, c, r]), esubd[c, r], γ_qidqim[c, r]))
            e_pia[c=comm, r =reg;δ_qia[c,r]], log(pia[c,r] .* qia[c,r]) == log(pid[c,r] * qid[c,r] + pim[c,r] * qim[c,r])
            e_psave[r=reg], log(psave[r]) == log(pinv[r] + sum(((qinv .- δ .* kb) .- qsave) .* pinv ./ sum(qinv .- δ .* kb)))

            # Trade - exports
            e_qms[c=comm, r=reg], log.(qms[c, r]) == log.(sum(qfm[c, :, r]) + qpm[c, r] + qgm[c, r] + qim[c, r])
            e_qxs[c=comm, r=reg], log.(Array(qxs[c, :, r])[δ_qxs[c, :, r]]) .== log.(ces(qms[c, r], Vector(pmds[c, :, r])[δ_qxs[c, :, r]], Vector(α_qxs[c, :, r])[δ_qxs[c, :, r]], esubm[c, r], γ_qxs[c, r]))
            e_pms[c=comm, r=reg], log(pms[c, r] * qms[c, r]) == log(sum(pmds[c, :, r] .* qxs[c, :, r]))

            # Trade - margins
            e_qtmfsd[m=marg, c=comm, s=reg, d=reg; δ_vtwr[m, c, s, d]], log(qtmfsd[m, c, s, d]) == log(α_qtmfsd[m, c, s, d] * qxs[c, s, d])
            e_ptrans[c=comm, s=reg, d=reg; δ_vtwr_sum[c, s, d]], log(ptrans[c, s, d] * sum(qtmfsd[:, c, s, d])) == log(sum(qtmfsd[:, c, s, d] .* pt[:]))
            e_qtm[m=marg], log(qtm[m]) == log(sum(qtmfsd[m, :, :, :]))
            e_qst[m=marg], log.(qst[m, :]) .== log.(ces(qtm[m], pds[m, :], Vector(α_qst[m, :]), esubs[m], γ_qst[m]))
            e_pt[m=marg], log.(pt[m] * qtm[m]) == log.(sum(pds[m, :] .* qst[m, :]))

            # Trade - imports / exports
            e_pfob[c=comm, s=reg, d=reg; δ_qxs[c, s, d]], log(pfob[c, s, d]) == log(pds[c, s] * tx[c, s] * txs[c, s, d])
            e_pcif[c=comm, s=reg, d=reg; δ_qxs[c, s, d]], log(pcif[c, s, d] * qxs[c, s, d]) == log(pfob[c, s, d] * (qxs[c, s, d]) + (δ_vtwr_sum[c, s, d] ? ptrans[c, s, d] * sum(Vector(qtmfsd[:, c, s, d])[δ_vtwr[:, c, s, d]]) : 0))
            e_pmds[c=comm, s=reg, d=reg; δ_qxs[c, s, d]], log(pmds[c, s, d]) == log(pcif[c, s, d] * tm[c, d] * tms[c, s, d])

            # Domestic market clearing 
            e_qds[c=comm, r=reg], log(qds[c, r]) == log(sum(qfd[c, :, r]) + qpd[c, r] + qgd[c, r] + qid[c, r])
            e_pds[c=comm, r=reg], log(qc[c, r]) == log(qds[c, r] + sum(qxs[c, r, :]) + (c ∈ marg ? qst[c, r] : 0))

            # Taxes
            e_pfd[c=comm, a=acts, r=reg], log(pfd[c, a, r]) == log(pds[c, r] * tfd[c, a, r])
            e_pfm[c=comm, a=acts, r=reg], log(pfm[c, a, r]) == log(pms[c, r] * tfm[c, a, r])
            e_ppd, log.(ppd) .== log.(pds .* tpd)
            e_ppm, log.(ppm) .== log.(pms .* tpm)
            e_pgd[c=comm, r=reg; δ_qga[c, r]], log(pgd[c, r]) == log(pds[c, r] * tgd[c, r])
            e_pgm[c=comm, r=reg; δ_qga[c, r]], log(pgm[c, r]) == log(pms[c, r] * tgm[c, r])
            e_pid[c=comm, r=reg; δ_qia[c, r]], log(pid[c, r]) == log(pds[c, r] * tid[c, r])
            e_pim[c=comm, r=reg; δ_qia[c, r]], log(pim[c, r]) == log(pms[c, r] * tim[c, r])

            # Factor Market
            e_pe1[e=endwm, r=reg], log.(qe[e, r]) == log.(sum(qfe[e, :, r]))
            e_qes1[e=endwm, a=acts, r=reg], log(pes[e, a, r]) == log(pe[e, r])
            e_qes2[e=endws, r=reg], log.(Vector(qes[e, :, r])[δ_evfp[e, :, r]]) .== log.(Vector(ces(qe[e, r], Vector(pes[e, :, r])[δ_evfp[e, :, r]], Vector(α_qes2[e, :, r])[δ_evfp[e, :, r]], etrae[e, r], γ_qes2[e, r])))
            e_pe2[e=endws, r=reg], log(pe[e, r] * qe[e, r]) == log(sum(pes[e, :, r] .* qes[e, :, r]))
            e_qes3[e=endwf, a=acts, r=reg; δ_evfp[e, a, r]], log(qes[e, a, r]) == log(qesf[e, a, r])
            e_pes[e=endw, a=acts, r=reg; δ_evfp[e, a, r]], log(peb[e, a, r]) == log(pes[e, a, r] * tinc[e, a, r])

            # Investment is a fixed share of global investment
            e_qinv, log.(qinv) .== log.(Vector(α_qinv) .* globalcgds .+ δ .* kb)
            e_pcgdswld, log(pcgdswld) == log(sum(pinv .* (qinv .- δ .* kb)) / sum(qinv .- δ .* kb))
            e_walras_sup, log(walras_sup) == log(pcgdswld * globalcgds)
            e_walras_dem, log(walras_dem) == log(sum(psave .* qsave))

            # Pfactwld
            e_rorg, log(pfactwld * sum(qfe)) == log(sum(peb .* qfe))

            # Capital accumulation
            e_kb[r=reg], log(ρ[r] * kb[r]) == log(sum(qe[endwc, r]))
            e_ke, log.(ke) .== log.(qinv .+ (1 .- δ) .* kb)

            # Values
            cvdfp, log.(vdfp) .== log.(pfd .* qfd)
            cvmfp, log.(vmfp) .== log.(pfm .* qfm)
            cvdpp, log.(vdpp) .== log.(ppd .* qpd)
            cvmpp, log.(vmpp) .== log.(ppm .* qpm)
            cvdgp[c=comm, r=reg; δ_qga[c, r]], log(vdgp[c, r]) == log(pgd[c, r] * qgd[c, r])
            cvmgp[c=comm, r=reg; δ_qga[c, r]], log(vmgp[c, r]) == log(pgm[c, r] * qgm[c, r])
            cvdip[c=comm, r=reg; δ_qia[c, r]], log(vdip[c, r]) .== log(pid[c, r] * qid[c, r])
            cvmip[c=comm, r=reg; δ_qia[c, r]], log(vmip[c, r]) .== log(pim[c, r] * qim[c, r])
            cevfp[e=endw, a=acts, r=reg; δ_evfp[e, a, r]], log.(evfp[e, a, r]) == log(pfe[e, a, r] * qfe[e, a, r])
            cevos[e=endw, a=acts, r=reg; δ_evfp[e, a, r]], log(evos[e, a, r]) == log(pes[e, a, r] * qes[e, a, r])
            cvfob[c=comm, s=reg, d=reg; δ_qxs[c, s, d]], log(vfob[c, s, d]) == log(pfob[c, s, d] * qxs[c, s, d])
            cvcif[c=comm, s=reg, d=reg; δ_qxs[c, s, d]], log(vcif[c, s, d]) == log(pcif[c, s, d] * qxs[c, s, d])
            cvst[m=marg], log.(vst[m, :]) .== log.(pds[m, :] .* qst[m, :])
            cvtwr[m=marg, c=comm, s=reg, d=reg; δ_vtwr[m, c, s, d]], log(vtwr[m, c, s, d]) == log.(pt[m] * qtmfsd[m, c, s, d])
            cmaks[c=comm, a=acts, r=reg; δ_maks[c, a, r]], log(maks[c, a, r]) == log(ps[c, a, r] * qca[c, a, r])
            cvkb, log.(vkb) .== log.(kb .* pinv)

            # Soft parameter constraints
            sf_α_qxs[c=comm, d=reg], log(sum(α_qxs[c, :, d])) == log(ϵ_qxs[c, d])
            sf_α_qfe[a=acts, r=reg], log(sum(α_qfe[:, a, r])) == log(ϵ_qfe[a, r])
            sf_α_qes2[e=endws, r=reg], log(sum(α_qes2[e, :, r])) == log(ϵ_qes2[e, r])
            sf_α_qfdqfm[c=comm, a=acts, r=reg], log(sum(α_qfdqfm[:, c, a, r])) == log(ϵ_qfdqfm[c, a, r])
            sf_α_qpdqpm[c=comm, r=reg], log(sum(α_qpdqpm[:, c, r])) == log(ϵ_qpdqpm[c, r])
            sf_α_qgdqgm[c=comm, r=reg; δ_qga[c, r]], log(sum(α_qgdqgm[:, c, r])) == log(ϵ_qgdqgm[c, r])
            sf_α_qidqim[c=comm, r=reg; δ_qia[c, r]], log(sum(α_qidqim[:, c, r])) == log(ϵ_qidqim[c, r])
            sf_α_qia[r=reg], log(sum(α_qia[:, r])) == log(ϵ_qia[r])
            sf_α_qga[r=reg], log(sum(α_qga[:, r])) == log(ϵ_qga[r])
            sf_α_qfa[a=acts, r=reg], log(sum(α_qfa[:, a, r])) == log(ϵ_qfa[a, r])
            sf_α_qintva[a=acts, r=reg], log(sum(α_qintva[["int", "va"], a, r])) == log(ϵ_qintva[a, r])
            sf_α_qinv, log(sum(α_qinv)) == log(ϵ_qinv)
            sf_save, log.(σsave .+ σyp .+ σyg) .== log(1)

            # Shares (helpers)
            e_σ_vp[c=comm, r=reg], log(σ_vp[c, r]) + log(sum(ppd[:, r] .* qpd[:, r] .+ ppm[:, r] .* qpm[:, r])) == log(ppd[c, r] .* qpd[c, r] + ppm[c, r] .* qpm[c, r])
            e_σ_vdp[c=comm, r=reg], log(σ_vdp[c, r]) + log(ppd[c, r] .* qpd[c, r] .+ ppm[c, r] .* qpm[c, r]) == log(ppd[c, r] .* qpd[c, r])
            e_σ_vg[c=comm, r=reg; δ_qga[c, r]], log(σ_vg[c, r]) + log(sum(pgd[:, r] .* qgd[:, r] .+ pgm[:, r] .* qgm[:, r])) == log(pgd[c, r] * qgd[c, r] + pgm[c, r] * qgm[c, r])
            e_σ_vdg[c=comm, r=reg; δ_qga[c, r]], log(σ_vdg[c, r]) + log(pgd[c, r] .* qgd[c, r] .+ pgm[c, r] .* qgm[c, r]) == log(pgd[c, r] .* qgd[c, r])
            e_σ_vi[c=comm, r=reg; δ_qia[c, r]], log(σ_vi[c, r]) + log(sum(pid[:, r] .* qid[:, r] .+ pim[:, r] .* qim[:, r])) == log(pid[c, r] .* qid[c, r] + pim[c, r] .* qim[c, r])
            e_σ_vdi[c=comm, r=reg; δ_qia[c, r]], log(σ_vdi[c, r]) + log(pid[c, r] .* qid[c, r] .+ pim[c, r] .* qim[c, r]) == log(pid[c, r] .* qid[c, r])
            e_σ_vf[c=comm, a=acts, r=reg], log(σ_vf[c, a, r]) + log(sum(pfd[:, a, r] .* qfd[:, a, r] .+ pfm[:, a, r] .* qfm[:, a, r])) == log(pfd[c, a, r] .* qfd[c, a, r] + pfm[c, a, r] .* qfm[c, a, r])
            e_σ_vdf[c=comm, a=acts, r=reg], log(σ_vdf[c, a, r]) + log(pfd[c, a, r] .* qfd[c, a, r] .+ pfm[c, a, r] .* qfm[c, a, r]) == log(pfd[c, a, r] .* qfd[c, a, r])
            e_σ_vif[a=acts, r=reg], log(σ_vif[a, r]) + log(pva[a, r] * qva[a, r] + pint[a, r] * qint[a, r]) == log(pint[a, r] * qint[a, r])
            e_σ_vff[e=endw, a=acts, r=reg; δ_evfp[e, a, r]], log(σ_vff[e, a, r]) + log(sum(pfe[:, a, r] .* qfe[:, a, r])) == log(pfe[e, a, r] * qfe[e, a, r])
            e_σ_vtwr[m=marg, c=comm, s=reg, d=reg; δ_vtwr[m, c, s, d]], log(σ_vtwr[m, c, s, d]) + log(pcif[c, s, d] * qxs[c, s, d]) == log(pt[m] * qtmfsd[m, c, s, d])
            e_σ_qxs[c=comm, s=reg, d=reg; δ_qxs[c, s, d]], log(σ_qxs[c, s, d]) + log(sum(pcif[c, :, d] .* qxs[c, :, d])) == log(pcif[c, s, d] .* qxs[c, s, d])
            e_σ_qinv[r=reg], σ_qinv[r] == psave[r] * qsave[r] / sum(psave .* qsave)
            e_σ_ρ[r=reg], log((kb[r] * pinv[r]) * σ_ρ[r]) == log(sum(pes[:, :, r] .* qes[:, :, r]))
        end
    )

    # Remove structurally missing variables
    # Remove qfe, qes, pfe, peb, α_qes2, α_qfe, cevos, cevfp if there is no use of the factor
    delete.(mc.model, Array(qfe)[.!δ_evfp])
    delete.(mc.model, Array(qes)[.!δ_evfp])
    delete.(mc.model, Array(pes)[.!δ_evfp])
    delete.(mc.model, Array(pfe)[.!δ_evfp])
    delete.(mc.model, Array(peb)[.!δ_evfp])
    delete.(mc.model, Array(α_qes2[endws, :, :])[.!δ_evfp[endws, :, :]])
    delete.(mc.model, Array(α_qfe)[.!δ_evfp])
    delete.(mc.model, Array(evos)[.!δ_evfp])
    delete.(mc.model, Array(evfp)[.!δ_evfp])
    delete.(mc.model, Array(σ_vff)[.!δ_evfp])

    # Remove maks when there is no output produced
    delete.(mc.model, Array(maks)[.!δ_maks])
    delete.(mc.model, Array(ps)[.!δ_maks])
    delete.(mc.model, Array(qca)[.!δ_maks])
    delete.(mc.model, Array(pca)[.!δ_maks])
    delete.(mc.model, Array(α_qca)[.!δ_maks])

    # Remove trade when no trade is allowed
    delete.(mc.model, Array(α_qxs)[.!δ_qxs])
    delete.(mc.model, Array(qxs)[.!δ_qxs])
    delete.(mc.model, Array(pmds)[.!δ_qxs])
    delete.(mc.model, Array(pcif)[.!δ_qxs])
    delete.(mc.model, Array(pfob)[.!δ_qxs])
    delete.(mc.model, Array(ptrans)[.!δ_qxs.||.!δ_vtwr_sum])
    delete.(mc.model, Array(vcif)[.!δ_qxs])
    delete.(mc.model, Array(vfob)[.!δ_qxs])
    delete.(mc.model, Array(σ_qxs)[.!δ_qxs])

    # Remove margins when no margins are allowed
    delete.(mc.model, Array(qtmfsd)[.!δ_vtwr])
    delete.(mc.model, Array(σ_vtwr)[.!δ_vtwr])
    delete.(mc.model, Array(α_qtmfsd)[.!δ_vtwr])
    delete.(mc.model, Array(vtwr)[.!δ_vtwr])

    # Remove government consumption if not present
    delete.(mc.model, Array(qga)[.!δ_qga])
    delete.(mc.model, Array(pga)[.!δ_qga])
    delete.(mc.model, Array(α_qga)[.!δ_qga])
    delete.(mc.model, Array(qgm)[.!δ_qga])
    delete.(mc.model, Array(pgm)[.!δ_qga])
    delete.(mc.model, Array(qgd)[.!δ_qga])
    delete.(mc.model, Array(pgd)[.!δ_qga])
    delete.(mc.model, Array(α_qgdqgm["dom",:,:])[.!δ_qga])
    delete.(mc.model, Array(α_qgdqgm["imp",:,:])[.!δ_qga])
    delete.(mc.model, Array(γ_qgdqgm[:,:])[.!δ_qga])
    delete.(mc.model, Array(σ_vdg)[.!δ_qga])
    delete.(mc.model, Array(σ_vg)[.!δ_qga])
    delete.(mc.model, Array(ϵ_qgdqgm)[.!δ_qga])
    delete.(mc.model, Array(tgm)[.!δ_qga])
    delete.(mc.model, Array(tgd)[.!δ_qga])
    delete.(mc.model, Array(vdgp)[.!δ_qga])
    delete.(mc.model, Array(vmgp)[.!δ_qga])
   
    # Remove investment consumption if not present
    delete.(mc.model, Array(qia)[.!δ_qia])
    delete.(mc.model, Array(pia)[.!δ_qia])
    delete.(mc.model, Array(α_qia)[.!δ_qia])
    delete.(mc.model, Array(qim)[.!δ_qia])
    delete.(mc.model, Array(pim)[.!δ_qia])
    delete.(mc.model, Array(qid)[.!δ_qia])
    delete.(mc.model, Array(pid)[.!δ_qia])
    delete.(mc.model, Array(α_qidqim["dom",:,:])[.!δ_qia])
    delete.(mc.model, Array(α_qidqim["imp",:,:])[.!δ_qia])
    delete.(mc.model, Array(γ_qidqim[:,:])[.!δ_qia])
    delete.(mc.model, Array(σ_vdi)[.!δ_qia])
    delete.(mc.model, Array(σ_vi)[.!δ_qia])
    delete.(mc.model, Array(ϵ_qidqim)[.!δ_qia])
    delete.(mc.model, Array(tim)[.!δ_qia])
    delete.(mc.model, Array(tid)[.!δ_qia])
    delete.(mc.model, Array(vdip)[.!δ_qia])
    delete.(mc.model, Array(vmip)[.!δ_qia])
   

    set_attribute(mc.model, "max_iter", max_iter)
    set_attribute(mc.model, "constr_viol_tol", constr_viol_tol)
    set_attribute(mc.model, "bound_push", bound_push)

    return nothing
end
