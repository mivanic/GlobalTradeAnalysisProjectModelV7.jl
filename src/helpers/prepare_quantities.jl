function prepare_quantities(; data, parameters, sets, hData)

    # Read all prices
    (; pint, pva, po, pfa, pfe, pfm, ps, pds, ppm, pes, pgm, pim, pcgdswld, psave, pmds, pt, pinv, walras_sup, ppa, ppd, pgd, pga, peb, tpd, pms, tpm, tgd, tgm, tid, tim, tfd, tfm, to, pfd, tfe, pfob, txs, pcif, tms, fincome, y, yg, yp) = NamedTuple(Dict(Symbol(k) => data[k] for k ∈ keys(data)))


    # Read all sets
    (; marg, comm, reg, endw, acts, endwf,endwm,endws) = NamedTuple(Dict(Symbol(k)=>sets[k] for k ∈ keys(sets)))

    vint = NamedArray(mapslices(sum, hData["vdfp"] .+ hData["vmfp"], dims=1)[1, :, :], names(hData["vdfp"])[[2, 3]])
    qint = vint ./ pint

    vva = NamedArray(mapslices(sum, hData["evfp"], dims=1)[1, :, :], names(hData["evfp"])[[2, 3]])
    qva = vva ./ pva

    vo = NamedArray(mapslices(sum, hData["maks"], dims=1)[1, :, :], names(hData["maks"])[[2, 3]])
    qo = vo ./ po

    vfa = hData["vdfp"] .+ hData["vmfp"]
    qfa = vfa ./ pfa

    vfe = hData["evfp"]
    qfe = vfe ./ pfe
    qfe[isnan.(qfe)].=0

    vfm = hData["vmfp"]
    qfm = vfm ./ pfm

    vfd = hData["vdfp"]
    qfd = vfd ./ pfd #repeat(reshape(pds, [length(comm), 1, length(reg)]...), inner=[1, length(acts), 1])

    vca = hData["maks"]
    qca = vca ./ ps

    vc = NamedArray(mapslices(sum, hData["makb"], dims=2)[:, 1, :], names(hData["maks"])[[1, 3]])
    qc = vc ./ pds

    vpd = hData["vdpp"]
    qpd = vpd ./ ppd

    vpm = hData["vmpp"]
    qpm = vpm ./ ppm

    qpa = (vpd .+ vpm) ./ ppa

    vgd = hData["vdgp"]
    qgd = vgd ./ pgd

    vgm = hData["vmgp"]
    qgm = vgm ./ pgm

    qga = (vgd .+ vgm) ./ pga

    vsave = hData["save"] #.+ hData["vdep"]
    qsave = vsave ./ psave

    vdip = hData["vdip"]
    vdib = hData["vdib"]
    qid = vdib ./ pds

    vim = hData["vmip"]
    qim = vim ./ pim

    qia = qid .+ qim

    vinv = NamedArray((mapslices(sum, vdip .+ vim, dims=1))[1, :], names(qid, 2))
    qinv = NamedArray(mapslices(sum, qia, dims=1)[1, :], reg) #vinv ./ pinv

    qms = NamedArray(mapslices(sum, hData["vmsb"], dims = [2])[:,1,:] ./ pms, (comm, reg)) #NamedArray(mapslices(sum, qfm, dims=2)[:, 1, :] .+ qpm .+ qgm .+ qim, names(qpm))

    vxs = hData["vmsb"]
    qxs = vxs ./ pmds

    vtmfsd = hData["vtwr"]
    qtmfsd = vtmfsd ./ repeat(reshape(pt, [length(marg), 1, 1, 1]...), inner=[1, length(comm), length(reg), length(reg)])

    qtm = NamedArray(mapslices(sum, qtmfsd, dims=[2, 3, 4])[:, 1, 1, 1], marg)

    vst = hData["vst"]
    qst = vst ./ pds[marg, :]

    qds = NamedArray(mapslices(sum, qfd, dims=2)[:, 1, :] .+ qpd .+ qgd .+ qid, names(qid))
    ves = hData["evos"]
    qes = ves ./ pes

    qes[isnan.(qes)] .= 0

    #qe = NamedArray(mapslices(sum, qes, dims=2)[:, 1, :], names(qes)[[1, 3]])[endwms, :]

    qesf = qes[endwf, acts, reg]

    globalcgds = sum(qinv)


    # Calculate values
    #yp = NamedArray(mapslices(sum, qpa .* ppa, dims=1)[1, :],reg)
    #yg = NamedArray(mapslices(sum, qga .* pga, dims=1)[1, :],reg)
    #fincome = NamedArray(mapslices(sum, peb .* qes, dims=[1, 2])[1, 1, :], reg)
    walras_sup = pcgdswld * globalcgds
    walras_dem = pcgdswld * globalcgds

    # y = fincome .+ mapslices(sum, qpd .* pds .* (tpd .- 1) .+ qpm .* pms .* (tpm .- 1) .+ qgd .* pds .* (tgd .- 1) .+ qgm .* pms .* (tgm .- 1) .+ qid .* pds .* (tid .- 1) .+ qim .* pms .* (tim .- 1), dims=1)[1, :] .+
    #     mapslices(sum, qfd .* pfd ./ tfd .* (tfd .- 1) .+ qfm .* pfm ./ tfm .* (tfm .- 1) .+ qca .* ps .* (to .- 1), dims=[1, 2])[1, 1, :] .+
    #     mapslices(sum, qfe .* peb .* (tfe .- 1), dims=[1, 2])[1, 1, :] .+
    #     mapslices(sum, qxs .* pfob ./ txs .* (txs .- 1), dims=[1, 3])[1, :, 1] .+
    #     mapslices(sum, qxs .* pcif .* (tms .- 1), dims=[1, 2])[1, 1, :]


    kb = hData["vkb"] ./ pinv
    ke = (hData["vkb"] .- hData["vdep"]) ./pinv .+ qinv 
   
    Φ = NamedArray(ones(length(reg)), reg)
    Φᴾ = NamedArray(ones(length(reg)), reg)

    quantities = Dict(
        :qint => qint,
        :qva => qva,
        :qo => qo,
        :qfa => qfa,
        :qfe => qfe,
        :qfm => qfm,
        :qfd => qfd,
        :qca => qca,
        :qc => qc,
        :qes => qfe,
        :qpd => qpd,
        :qpm => qpm,
        :qpa => qpa,
        :qgd => qgd,
        :qgm => qgm,
        :qga => qga,
        :qsave => qsave,
        :qinv => qinv,
        :qid => qid,
        :qim => qim,
        :qia => qia,
        :qms => qms,
        :qxs => qxs,
        :qtmfsd => qtmfsd,
        :qtm => qtm,
        :qst => qst,
        :qds => qds,
        :qes => qes,
        #:qe => qe,
        :globalcgds => globalcgds,
        :qesf => qesf,
        :yp => yp,
        :yg => yg,
        :y => y,
        :fincome => fincome,
        :walras_sup => walras_sup,
        :walras_dem => walras_dem,
        :kb => kb,
        :ke => ke,
        :evfp => hData["evfp"],
        :maks => hData["maks"],
        :vtwr => hData["vtwr"],
        :vcif => hData["vcif"],
        :evos => hData["evos"],
        :vdfp => hData["vdfp"],
        :vmfp => hData["vmfp"],
        :vdpp => hData["vdpp"],
        :vmpp => hData["vmpp"],
        :vdgp => hData["vdgp"],
        :vmgp => hData["vmgp"],
        :vdip => hData["vdip"],
        :vmip => hData["vmip"],
        :Φᴾ => Φᴾ,
        :Φ => Φ
    )

    return merge(data, Dict(String(k)=>quantities[k] for k ∈ keys(quantities)))
end