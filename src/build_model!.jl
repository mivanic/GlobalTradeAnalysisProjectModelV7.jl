function build_model!(mc; max_iter=50, constr_viol_tol=1e-8, bound_push=1e-15)

    # Structural parameters (some CES/CET options are not happening)
    δ_evfp = .!isnan.(mc.data["α_qfe"]) .&& mc.data["α_qfe"] .!= 0
    δ_maks = .!isnan.(mc.data["α_qca"]) .&& mc.data["α_qca"] .!= 0
    δ_vtwr = .!isnan.(mc.data["α_qtmfsd"]) .&& mc.data["α_qtmfsd"] .!= 0
    δ_vtwr_sum = NamedArray(.!mapslices(all, isnan.(mc.data["α_qtmfsd"]) .|| mc.data["α_qtmfsd"] .== 0, dims=1)[1, :, :, :], names(mc.data["vtwr"])[[2, 3, 4]])
    δ_qxs = .!isnan.(mc.data["α_qxs"]) .&& mc.data["α_qxs"] .!= 0

    # Read  sets
    (; reg, comm, marg, acts, endw, endwc, endws, endwm, endwms, endwf) = NamedTuple(Dict(Symbol(k) => mc.sets[k] for k ∈ keys(mc.sets)))

    # Read hard parameters
    (; esubt, esubc, esubva, esubd, etraq, esubq, subpar, incpar, etrae, esubg, esubm, esubs) = NamedTuple(Dict(Symbol(k) => mc.parameters[k] for k ∈ keys(mc.parameters)))

    # Set up the general constraints
    p_min = 1e-8
    p_max = 1e+8
    q_min = 1e-8
    q_max = 1e+12
    y_min = 1e-8
    y_max = 1e+12
    t_min = -1e2
    t_max = 1e2

    # All variables used in the model
    @variables(mc.model,
        begin
            # Population
            q_min <= pop[reg] <= q_max

            # Firms top nest
            q_min <= qint[acts, reg] <= q_max
            p_min <= pint[acts, reg] <= p_max
            q_min <= qva[acts, reg] <= q_max
            p_min <= pva[acts, reg] <= p_max
            q_min <= qo[acts, reg] <= q_max
            p_min <= po[acts, reg] <= p_max

            # Firms second nest
            q_min <= qfa[comm, acts, reg] <= q_max
            p_min <= pfa[comm, acts, reg] <= p_max
            q_min <= qfe[endw, acts, reg] <= q_max
            p_min <= pfe[endw, acts, reg] <= p_max
            p_min <= tfe[endw, acts, reg] <= p_max
            q_min <= qfd[comm, acts, reg] <= q_max
            p_min <= pfd[comm, acts, reg] <= p_max
            q_min <= qfm[comm, acts, reg] <= q_max
            p_min <= pfm[comm, acts, reg] <= p_max

            # # # Firm distribution
            q_min <= qca[comm, acts, reg] <= q_max
            p_min <= pca[comm, acts, reg] <= p_max
            p_min <= ps[comm, acts, reg] <= p_max
            q_min <= qc[comm, reg] <= q_max
            p_min <= pds[comm, reg] <= p_max
            t_min <= to[comm, acts, reg] <= t_max

            # Endowments
            p_min <= peb[endw, acts, reg] <= p_max
            q_min <= qes[endw, acts, reg] <= q_max
            p_min <= pfactor[reg] <= p_max

            # Income
            y_min <= fincome[reg] <= y_max
            y_min <= y[reg] <= y_max

            q_min <= u[reg] <= q_max
            q_min <= ug[reg] <= q_max
            q_min <= us[reg] <= q_max

            # Private consumption        
            y_min <= yp[reg] <= y_max
            q_min <= up[reg] <= q_max
            0 <= uelas[reg] <= 10
            0 <= uepriv[reg] <= 10
            p_min <= ppa[comm, reg] <= p_max
            q_min <= qpa[comm, reg] <= q_max
            p_min <= ppd[comm, reg] <= p_max
            q_min <= qpd[comm, reg] <= q_max
            p_min <= ppm[comm, reg] <= p_max
            q_min <= qpm[comm, reg] <= q_max
            p_min <= ppriv[reg] <= p_max

            # Government consumption
            y_min <= yg[reg] <= y_max
            p_min <= pgov[reg] <= p_max
            p_min <= pga[comm, reg] <= p_max
            q_min <= qga[comm, reg] <= q_max
            p_min <= pgd[comm, reg] <= p_max
            q_min <= qgd[comm, reg] <= q_max
            p_min <= pgm[comm, reg] <= p_max
            q_min <= qgm[comm, reg] <= q_max

            # Saving
            -q_max <= qsave[reg] <= q_max
            p_min <= psave[reg] <= p_max

            # Investment consumption
            p_min <= pia[comm, reg] <= p_max
            q_min <= qia[comm, reg] <= q_max
            p_min <= pid[comm, reg] <= p_max
            q_min <= qid[comm, reg] <= q_max
            p_min <= pim[comm, reg] <= p_max
            q_min <= qim[comm, reg] <= q_max
            q_min <= qinv[reg] <= q_max
            p_min <= pinv[reg] <= p_max

            # Trade - exports
            q_min <= qms[comm, reg] <= q_max
            q_min <= qxs[comm, reg, reg] <= q_max
            p_min <= pmds[comm, reg, reg] <= p_max
            p_min <= pms[comm, reg] <= p_max

            # Trade - margins
            p_min <= ptrans[comm, reg, reg] <= p_max
            q_min <= qtmfsd[marg, comm, reg, reg] <= q_max
            p_min <= pt[marg] <= p_max
            q_min <= qtm[marg] <= q_max
            q_min <= qst[marg, reg] <= q_max

            # Trade - imports/exports
            t_min <= txs[comm, reg, reg] <= t_max
            t_min <= tx[comm, reg] <= t_max
            p_min <= pfob[comm, reg, reg] <= p_max
            p_min <= pcif[comm, reg, reg] <= p_max
            t_min <= tms[comm, reg, reg] <= t_max
            t_min <= tm[comm, reg] <= t_max

            # Domestic market clearing
            q_min <= qds[comm, reg] <= q_max

            # Taxes
            t_min <= tfd[comm, acts, reg] <= t_max
            t_min <= tfm[comm, acts, reg] <= t_max
            t_min <= tpd[comm, reg] <= t_max
            t_min <= tpm[comm, reg] <= t_max
            t_min <= tgd[comm, reg] <= t_max
            t_min <= tgm[comm, reg] <= t_max
            t_min <= tid[comm, reg] <= t_max
            t_min <= tim[comm, reg] <= t_max

            #Factor Market
            p_min <= pes[endw, acts, reg] <= p_max
            p_min <= pe[endwms, reg] <= p_max
            q_min <= qe[endwms, reg] <= q_max
            q_min <= qesf[endwf, acts, reg] <= q_max
            p_min <= tinc[endw, acts, reg] <= p_max

            # Global Investment
            1 <= globalcgds <= q_max
            p_min <= pcgdswld <= p_max
            q_min <= walras_sup <= q_max
            q_min <= walras_dem <= q_max
            p_min <= pfactwld <= p_max

            # Capital stocks
            q_min <= kb[reg] <= q_max
            q_min <= ke[reg] <= q_max

            # Soft parameters
            1e-8 <= α_qintva[["int", "va"], acts, reg] <= 1
            1e-8 <= ϵ_qintva[acts, reg]
            1e-8 <= γ_qintva[acts, reg]
            1e-8 <= α_qfa[comm, acts, reg] <= 1
            1e-8 <= ϵ_qfa[acts, reg]
            1e-8 <= γ_qfa[acts, reg]
            1e-8 <= α_qfe[endw, acts, reg] <= 1
            1e-8 <= γ_qfe[acts, reg]
            0 <= ϵ_qfe[acts, reg]
            1e-8 <= α_qfdqfm[["dom", "imp"], comm, acts, reg] <= 1
            1e-8 <= γ_qfdqfm[comm, acts, reg]
            0 <= ϵ_qfdqfm[comm, acts, reg]
            1e-8 <= α_qca[comm, acts, reg] <= 1
            1e-8 <= γ_qca[acts, reg]
            1e-8 <= α_pca[comm, acts, reg] <= 1
            1e-8 <= γ_pca[acts, reg]
            1e-8 <= σyp[reg] <= 1
            1e-8 <= σyg[reg] <= 1
            -1 <= σsave[reg] <= 1
            1e-8 <= α_qga[comm, reg] <= 1
            1e-8 <= ϵ_qga[reg]
            1e-8 <= γ_qga[reg]
            1e-8 <= α_qia[comm, reg] <= 1
            1e-8 <= γ_qia[reg]
            1e-8 <= ϵ_qia[reg]
            1e-8 <= α_qpdqpm[["dom", "imp"], acts, reg] <= 1
            1e-8 <= γ_qpdqpm[acts, reg]
            0 <= ϵ_qpdqpm[acts, reg]
            1e-8 <= α_qgdqgm[["dom", "imp"], acts, reg] <= 1
            1e-8 <= γ_qgdqgm[acts, reg]
            0 <= ϵ_qgdqgm[acts, reg]
            1e-8 <= β_qpa[comm, reg]
            1e-8 <= α_qidqim[["dom", "imp"], acts, reg] <= 1
            1e-8 <= γ_qidqim[acts, reg]
            0 <= ϵ_qidqim[acts, reg]
            1e-8 <= α_qxs[comm, reg, reg] <= 1
            1e-8 <= γ_qxs[comm, reg]
            0 <= ϵ_qxs[comm, reg]
            1e-8 <= α_qtmfsd[marg, comm, reg, reg] <= 1
            1e-8 <= α_qst[marg, reg] <= 1
            1e-8 <= γ_qst[marg]
            1e-8 <= α_qes2[endws, acts, reg] <= 1
            1e-8 <= γ_qes2[endws, reg]
            1e-8 <= ϵ_qes2[endws, reg]
            -1 <= α_qinv[reg] <= 1
            0 <= ϵ_qinv

            0 <= δ[reg] <= 1
            0 <= ρ[reg] <= 1

            # Values
            0 <= vdfp[comm, acts, reg]
            0 <= vmfp[comm, acts, reg]
            0 <= vdpp[comm, reg]
            0 <= vmpp[comm, reg]
            0 <= vdgp[comm, reg]
            0 <= vmgp[comm, reg]
            0 <= vdip[comm, reg]
            0 <= vmip[comm, reg]
            0 <= evfp[endw, acts, reg]
            0 <= evos[endw, acts, reg]
            0 <= vfob[comm, reg, reg]
            0 <= vcif[comm, reg, reg]
            0 <= vst[marg, reg]
            0 <= vtwr[marg, comm, reg, reg]
            0 <= maks[comm, acts, reg]
            0 <= vkb[reg]

            # Shares (helpers to calibrate)
            0 <= σ_vp[comm, reg]
            0 <= σ_vdp[comm, reg]
            0 <= σ_vg[comm, reg]
            0 <= σ_vdg[comm, reg]
            0 <= σ_vi[comm, reg]
            0 <= σ_vdi[comm, reg]
            0 <= σ_vf[comm, acts, reg]
            0 <= σ_vdf[comm, acts, reg]
            0 <= σ_vff[endw, acts, reg]
            0 <= σ_vif[acts, reg]
            0 <= σ_vtwr[marg, comm, reg, reg]
            0 <= σ_qxs[comm, reg, reg]
            -1 <= σ_qinv[reg] <= 1
            0 <= σ_ρ[reg]
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

    # All model equations
    @constraints(mc.model,
        begin
            # Firms (top nest)
            e_qintva[a=acts, r=reg], log.([qint[a, r], qva[a, r]]) .== log.(ces(qo[a, r], [pint[a, r], pva[a, r]], α_qintva[:, a, r], esubt[a, r], γ_qintva[a, r]))
            e_qo, log.(qo .* po) .== log.(qva .* pva .+ qint .* pint)

            # Firms (second nest)
            e_qfa[a=acts, r=reg], log.(qfa[:, a, r]) .== log.(ces(qint[a, r], pfa[:, a, r], Vector(α_qfa[:, a, r]), esubc[a, r], γ_qfa[a, r]))
            e_pint[a=acts, r=reg], log.(qint[a, r] * pint[a, r]) == log.(sum(pfa[:, a, r] .* qfa[:, a, r]))
            e_qfe[a=acts, r=reg], log.(Vector(qfe[:, a, r])[δ_evfp[:, a, r]]) .== log.(Vector(ces(qva[a, r], Vector(pfe[:, a, r])[δ_evfp[:, a, r].!=0], Vector(α_qfe[:, a, r])[δ_evfp[:, a, r].!=0], esubva[a, r], γ_qfe[a, r])))
            e_pva[a=acts, r=reg], log.(qva[a, r] * pva[a, r]) == log.(sum(Vector(pfe[:, a, r] .* qfe[:, a, r])[δ_evfp[:, a, r].!=0]))
            e_qfdqfm[c=comm, a=acts, r=reg], log.([qfd[c, a, r], qfm[c, a, r]]) .== log.(ces(qfa[c, a, r], [pfd[c, a, r], pfm[c, a, r]], α_qfdqfm[:, c, a, r], esubd[c, r], γ_qfdqfm[c, a, r]))
            e_pfa, log.(pfa .* qfa) .== log.(qfd .* pfd .+ qfm .* pfm)

            # Firms (distribution)
            e_qca[a=acts, r=reg], log.(Vector(qca[:, a, r])[δ_maks[:, a, r]]) .== log.(Vector(ces(qo[a, r], Vector(ps[:, a, r])[δ_maks[:, a, r]], Vector(α_qca[:, a, r])[δ_maks[:, a, r]], etraq[a, r], γ_qca[a, r])))
            e_po[a=acts, r=reg], log.(po[a, r] * qo[a, r]) == log.(sum(Vector(qca[:, a, r] .* ps[:, a, r])[δ_maks[:, a, r]]))
            e_pca[c=comm, r=reg], log.((esubq[c, r] == 0 ? Vector(pca[c, :, r])[δ_maks[c, :, r]] : Vector(qca[c, :, r])[δ_maks[c, :, r]])) .== log.((esubq[c, r] == 0 ? pds[c, r] : Vector(ces(qc[c, r], Vector(pca[c, :, r])[δ_maks[c, :, r]], Vector(α_pca[c, :, r])[δ_maks[c, :, r]], 1 / esubq[c, r], γ_pca[c, r]))))
            e_qc[c=comm, r=reg], log.(pds[c, r] * qc[c, r]) == log.(sum(Array(pca[c, :, r] .* qca[c, :, r])[δ_maks[c, :, r]]))
            e_ps, log.(pca) .== log.(ps .* to)

            # Endowments
            e_peb[e=endw, a=acts, r=reg], log.(qfe[e, a, r]) == log.(qes[e, a, r])
            e_pfe[e=endw, a=acts, r=reg], log.(pfe[e, a, r]) == log.(peb[e, a, r] .* tfe[e, a, r])
            e_pfactor[r=reg], log(pfactor[r] * sum(Array(qfe[:, :, r])[δ_evfp[:, :, r]])) == log(sum(Array(qfe[:, :, r] .* peb[:, :, r])[δ_evfp[:, :, r]]))

            # Income
            e_fincome[r=reg], log.(fincome[r]) == log.(sum(Array(peb[:, :, r] .* qes[:, :, r])[δ_evfp[:, :, r]]) .- δ[r] .* pinv[r] .* kb[r])
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
                sum(Array(qca[:, :, r] .* ps[:, :, r] .* (to[:, :, r] .- 1))[δ_maks[:, :, r]]) +
                sum(Array(qfe[:, :, r] .* peb[:, :, r] .* (tfe[:, :, r] .- 1))[δ_evfp[:, :, r]]) +
                sum(Array(qxs[:, r, :] .* pfob[:, r, :] ./ txs[:, r, :] .* (txs[:, r, :] .- 1))[δ_qxs[:, r, :]]) +
                sum(Array(qxs[:, :, r] .* pcif[:, :, r] .* (tms[:, :, r] .- 1))[δ_qxs[:, :, r]])
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
            e_yg, log.(yg) .== log.(y .* Vector(σyg) .* uelas)

            # Government expenditure
            e_qga[r=reg], log.(pga[:, r] .* qga[:, r]) .== log.(yg[r] .* Vector(α_qga[:, r])) ##This one
            e_pgov[r=reg], log.(pgov[r] * sum(qga[:, r])) == log.(sum(qga[:, r] .* pga[:, r]))
            e_qgdqgm[c=comm, r=reg], log.([qgd[c, r], qgm[c, r]]) .== log.(ces(qga[c, r], [pgd[c, r], pgm[c, r]], α_qgdqgm[:, c, r], esubd[c, r], γ_qgdqgm[c, r]))
            e_pga, log.(qga .* pga) .== log.(pgd .* qgd .+ pgm .* qgm)

            # Saving
            e_qsave, log.(y .* uelas) .== log.(σyp .* y .* uelas .+ σyg .* y .* uelas .+ psave .* qsave)

            # Investment consumption
            e_qia[r=reg], log.(qia[:, r]) .== log.(ces(qinv[r], pia[:, r], Vector(α_qia[:, r]), 0, γ_qia[r]))
            e_pinv[r=reg], log.(pinv[r] * sum(qia[:, r])) == log.(sum(pia[:, r] .* qia[:, r]))
            e_qidqim[c=comm, r=reg], log.([qid[c, r], qim[c, r]]) .== log.(ces(qia[c, r], [pid[c, r], pim[c, r]], Vector(α_qidqim[:, c, r]), esubd[c, r], γ_qidqim[c, r]))
            e_pia, log.(pia .* qia) .== log.(pid .* qid .+ pim .* qim)
            #e_psave, log.(psave) .== log.(pinv)
            e_psave[r=reg], log(psave[r]) == log(pinv[r] + sum(((qinv .- δ .* kb) .- qsave) .* pinv ./ sum(qinv .- δ .* kb)))

            # Trade - exports
            e_qms[c=comm, r=reg], log.(qms[c, r]) == log.(sum(qfm[c, :, r]) + qpm[c, r] + qgm[c, r] + qim[c, r])
            e_qxs[c=comm, r=reg], log.(Array(qxs[c, :, r])[δ_qxs[c, :, r]]) .== log.(ces(qms[c, r], Vector(pmds[c, :, r])[δ_qxs[c, :, r]], Vector(α_qxs[c, :, r])[δ_qxs[c, :, r]], esubm[c, r], γ_qxs[c, r]))
            e_pms[c=comm, r=reg], log.(pms[c, r] * qms[c, r]) == log.(sum(Vector(pmds[c, :, r] .* qxs[c, :, r])[δ_qxs[c, :, r]]))

            # Trade - margins
            e_qtmfsd[m=marg, c=comm, s=reg, d=reg], log.(qtmfsd[m, c, s, d]) .== log.(α_qtmfsd[m, c, s, d] .* qxs[c, s, d])
            e_ptrans[c=comm, s=reg, d=reg], log.(ptrans[c, s, d] * sum(Vector(qtmfsd[:, c, s, d])[δ_vtwr[:, c, s, d]]; init=0.0)) == log.(sum(Vector(qtmfsd[:, c, s, d])[δ_vtwr[:, c, s, d]] .* Vector(pt[:]); init=0.0))
            e_qtm[m=marg], log.(qtm[m]) == log.(sum(Array(qtmfsd[m, :, :, :])[δ_vtwr[m, :, :, :]]))
            e_qst[m=marg], log.(qst[m, :]) .== log.(ces(qtm[m], pds[m, :], Vector(α_qst[m, :]), esubs[m], γ_qst[m]))
            e_pt[m=marg], log.(pt[m] * qtm[m]) == log.(sum(pds[m, :] .* qst[m, :]))

            # Trade - imports / exports
            e_pfob[c=comm, s=reg, d=reg], log(pfob[c, s, d]) == log(pds[c, s] * tx[c, s] * txs[c, s, d])
            e_pcif[c=comm, s=reg, d=reg], log(pcif[c, s, d] * qxs[c, s, d]) == log(pfob[c, s, d] * (qxs[c, s, d]) + (δ_vtwr_sum[c, s, d] ? ptrans[c, s, d] * sum(Vector(qtmfsd[:, c, s, d])[δ_vtwr[:, c, s, d]]) : 0))
            e_pmds[c=comm, s=reg, d=reg], log(pmds[c, s, d]) == log(pcif[c, s, d] * tm[c, d] * tms[c, s, d])

            # Domestic market clearing 
            e_qds[c=comm, r=reg], log(qds[c, r]) == log(sum(qfd[c, :, r]) + qpd[c, r] + qgd[c, r] + qid[c, r])
            e_pds[c=comm, r=reg], log(qc[c, r]) == log(qds[c, r] + sum(Vector(qxs[c, r, :])[δ_qxs[c, r, :]]) + (c ∈ marg ? qst[c, r] : 0))

            # Taxes
            e_pfd[c=comm, a=acts, r=reg], log(pfd[c, a, r]) == log(pds[c, r] * tfd[c, a, r])
            e_pfm[c=comm, a=acts, r=reg], log(pfm[c, a, r]) == log(pms[c, r] * tfm[c, a, r])
            e_ppd, log.(ppd) .== log.(pds .* tpd)
            e_ppm, log.(ppm) .== log.(pms .* tpm)
            e_pgd, log.(pgd) .== log.(pds .* tgd)
            e_pgm, log.(pgm) .== log.(pms .* tgm)
            e_pid, log.(pid) .== log.(pds .* tid)
            e_pim, log.(pim) .== log.(pms .* tim)

            # Factor Market
            e_pe1[e=endwm, r=reg], log.(qe[e, r]) == log.(sum(Vector(qfe[e, :, r])[δ_evfp[e, :, r]]))
            e_qes1[e=endwm, a=acts, r=reg], log(pes[e, a, r]) == log(pe[e, r])
            e_qes2[e=endws, r=reg], log.(Vector(qes[e, :, r])[δ_evfp[e, :, r]]) .== log.(Vector(ces(qe[e, r], Vector(pes[e, :, r])[δ_evfp[e, :, r]], Vector(α_qes2[e, :, r])[δ_evfp[e, :, r]], etrae[e, r], γ_qes2[e, r])))
            e_pe2[e=endws, r=reg], log(pe[e, r] * qe[e, r]) == log(sum(Vector(pes[e, :, r] .* qes[e, :, r])[δ_evfp[e, :, r]]))
            e_qes3[e=endwf, a=acts, r=reg], log(qes[e, a, r]) == log(qesf[e, a, r])
            e_pes[e=endw, a=acts, r=reg], log(peb[e, a, r]) == log(pes[e, a, r] * tinc[e, a, r])

            # Investment is a fixed share of global investment
            e_qinv, log.(qinv) .== log.(Vector(α_qinv) .* globalcgds .+ δ .* kb)
            e_pcgdswld, log(pcgdswld) == log(sum(pinv .* (qinv .- δ .* kb)) / sum(qinv .- δ .* kb))
            e_walras_sup, log(walras_sup) == log(pcgdswld * globalcgds)
            e_walras_dem, log(walras_dem) == log(sum(psave .* qsave))

            # Pfactwld
            e_rorg, log(pfactwld * sum(Array(qfe)[δ_evfp])) == log(sum(Array(peb .* qfe)[δ_evfp]))

            # Capital accumulation
            #e_kb[r=reg], log(ρ[r] * pinv[r] * kb[r]) == log(sum(qe[endwc, r] .* pe[endwc, r]))
            e_kb[r=reg], log(ρ[r] * kb[r]) == log(sum(qe[endwc, r]))
            e_ke, log.(ke) .== log.(qinv .+ (1 .- δ) .* kb)

            # Values
            cvdfp, log.(vdfp) .== log.(pfd .* qfd)
            cvmfp, log.(vmfp) .== log.(pfm .* qfm)
            cvdpp, log.(vdpp) .== log.(ppd .* qpd)
            cvmpp, log.(vmpp) .== log.(ppm .* qpm)
            cvdgp, log.(vdgp) .== log.(pgd .* qgd)
            cvmgp, log.(vmgp) .== log.(pgm .* qgm)
            cvdip, log.(vdip) .== log.(pid .* qid)
            cvmip, log.(vmip) .== log.(pim .* qim)
            cevfp, log.(Array(evfp)) .== log.(Array(pfe .* qfe))
            cevos, log.(Array(evos)) .== log.(Array(pes .* qes))
            cvfob, log.(vfob) .== log.(pfob .* qxs)
            cvcif, log.(vcif) .== log.(pcif .* qxs)
            cvst[m=marg], log.(vst[m, :]) .== log.(pds[m, :] .* qst[m, :])
            cvtwr[c=comm, s=reg, d=reg], log.(Vector(vtwr[:, c, s, d])[δ_vtwr[:, c, s, d]]) .== log.(Vector(pt .* qtmfsd[:, c, s, d])[δ_vtwr[:, c, s, d]])
            cmaks, log.(maks) .== log.(ps .* qca)
            cvkb, log.(vkb) .== log.(kb .* pinv)

            # Soft parameter constraints
            sf_α_qxs[c=comm, d=reg], log(sum(Vector(α_qxs[c, :, d])[δ_qxs[c, :, d]])) == log(ϵ_qxs[c, d])
            sf_α_qfe[a=acts, r=reg], log(sum(Vector(α_qfe[:, a, r])[δ_evfp[:, a, r]])) == log(ϵ_qfe[a, r])
            sf_α_qes2[e=endws, r=reg], log(sum(Vector(α_qes2[e, :, r])[δ_evfp[e, :, r]])) == log(ϵ_qes2[e, r])
            sf_α_qfdqfm[c=comm, a=acts, r=reg], log(sum(α_qfdqfm[:, c, a, r])) == log(ϵ_qfdqfm[c, a, r])
            sf_α_qpdqpm[c=comm, r=reg], log(sum(α_qpdqpm[:, c, r])) == log(ϵ_qpdqpm[c, r])
            sf_α_qgdqgm[c=comm, r=reg], log(sum(α_qgdqgm[:, c, r])) == log(ϵ_qgdqgm[c, r])
            sf_α_qidqim[c=comm, r=reg], log(sum(α_qidqim[:, c, r])) == log(ϵ_qidqim[c, r])
            sf_α_qia[r=reg], log(sum(α_qia[:, r])) == log(ϵ_qia[r])
            sf_α_qga[r=reg], log(sum(α_qga[:, r])) == log(ϵ_qga[r])
            sf_α_qfa[a=acts, r=reg], log(sum(α_qfa[:, a, r])) == log(ϵ_qfa[a, r])
            sf_α_qintva[a=acts, r=reg], log(sum(α_qintva[["int", "va"], a, r])) == log(ϵ_qintva[a, r])
            sf_α_qinv, log(sum(α_qinv)) == log(ϵ_qinv)
            sf_save, log.(σsave .+ σyp .+ σyg) .== log(1)

            # Shares (helpers)
            e_σ_vp[c=comm, r=reg], σ_vp[c, r] * sum(ppd[:, r] .* qpd[:, r] .+ ppm[:, r] .* qpm[:, r]) == ppd[c, r] .* qpd[c, r] + ppm[c, r] .* qpm[c, r]
            e_σ_vdp[c=comm, r=reg], σ_vdp[c, r] * (ppd[c, r] .* qpd[c, r] .+ ppm[c, r] .* qpm[c, r]) == ppd[c, r] .* qpd[c, r]
            e_σ_vg[c=comm, r=reg], σ_vg[c, r] * sum(pgd[:, r] .* qgd[:, r] .+ pgm[:, r] .* qgm[:, r]) == pgd[c, r] .* qgd[c, r] + pgm[c, r] .* qgm[c, r]
            e_σ_vdg[c=comm, r=reg], σ_vdg[c, r] * (pgd[c, r] .* qgd[c, r] .+ pgm[c, r] .* qgm[c, r]) == pgd[c, r] .* qgd[c, r]
            e_σ_vi[c=comm, r=reg], σ_vi[c, r] * sum(pid[:, r] .* qid[:, r] .+ pim[:, r] .* qim[:, r]) == pid[c, r] .* qid[c, r] + pim[c, r] .* qim[c, r]
            e_σ_vdi[c=comm, r=reg], σ_vdi[c, r] * (pid[c, r] .* qid[c, r] .+ pim[c, r] .* qim[c, r]) == pid[c, r] .* qid[c, r]
            e_σ_vf[c=comm, a=acts, r=reg], σ_vf[c, a, r] * sum(pfd[:, a, r] .* qfd[:, a, r] .+ pfm[:, a, r] .* qfm[:, a, r]) == pfd[c, a, r] .* qfd[c, a, r] + pfm[c, a, r] .* qfm[c, a, r]
            e_σ_vdf[c=comm, a=acts, r=reg], σ_vdf[c, a, r] * (pfd[c, a, r] .* qfd[c, a, r] .+ pfm[c, a, r] .* qfm[c, a, r]) == pfd[c, a, r] .* qfd[c, a, r]
            e_σ_vif[a=acts, r=reg], σ_vif[a, r] * (pva[a, r] * qva[a, r] + pint[a, r] * qint[a, r]) == pint[a, r] * qint[a, r]
            e_σ_vff[e=endw, a=acts, r=reg], σ_vff[e, a, r] * sum(Vector(pfe[:, a, r] .* qfe[:, a, r])[δ_evfp[:, a, r]]) == pfe[e, a, r] .* qfe[e, a, r]
            e_σ_vtwr[m=marg, c=comm, s=reg, d=reg], σ_vtwr[m, c, s, d] * pcif[c, s, d] * qxs[c, s, d] == pt[m] * qtmfsd[m, c, s, d]
            e_σ_qxs[c=comm, s=reg, d=reg], σ_qxs[c, s, d] * sum(Vector(pcif[c, :, d] .* qxs[c, :, d])[δ_qxs[c, :, d]]) == pcif[c, s, d] .* qxs[c, s, d]
            e_σ_qinv[r=reg], σ_qinv[r] == psave[r] * qsave[r] / sum(psave .* qsave)
            e_σ_ρ[r=reg], log((kb[r] * pinv[r]) * σ_ρ[r]) == log(sum(Array(pes[:, :, r] .* qes[:, :, r])[δ_evfp[:, :, r]]))
        end
    )

    # Remove equations on structural grounds

    # Factors not used
    delete.(mc.model, Array(cevos)[.!δ_evfp])
    delete.(mc.model, Array(cevfp)[.!δ_evfp])
    delete.(mc.model, Array(e_pfe)[.!δ_evfp])
    delete.(mc.model, Array(e_peb)[.!δ_evfp])
    delete.(mc.model, Array(e_pes)[.!δ_evfp])
    delete.(mc.model, Array(e_qes1)[.!δ_evfp[endwm, :, :]])
    delete.(mc.model, Array(e_qes3)[.!δ_evfp[endwf, :, :]])
    delete.(mc.model, Array(e_σ_vff)[.!δ_evfp])

    # Outputs not produced
    delete.(mc.model, Array(cmaks)[.!δ_maks])
    delete.(mc.model, Array(e_ps)[.!δ_maks])

    # Trade not allowed
    delete.(mc.model, Array(e_pmds)[.!δ_qxs])
    delete.(mc.model, Array(e_pcif)[.!δ_qxs])
    delete.(mc.model, Array(e_pfob)[.!δ_qxs])
    delete.(mc.model, Array(e_ptrans)[.!δ_qxs])
    delete.(mc.model, Array(cvcif)[.!δ_qxs])
    delete.(mc.model, Array(cvfob)[.!δ_qxs])
    delete.(mc.model, Array(e_σ_qxs)[.!δ_qxs])


    # Margins not allowed

    delete.(mc.model, Array(e_qtmfsd)[.!δ_vtwr])
    delete.(mc.model, Array(e_σ_vtwr)[.!δ_vtwr])
    delete.(mc.model, Array(e_ptrans)[δ_qxs.&&.!δ_vtwr_sum])


    set_attribute(mc.model, "max_iter", max_iter)
    set_attribute(mc.model, "constr_viol_tol", constr_viol_tol)
    set_attribute(mc.model, "bound_push", bound_push)

    return nothing
end
