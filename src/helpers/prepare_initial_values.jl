# This calculates some initial values 
function prepare_initial_values(; sets, hData, hParameters)
    # Read the sets
    (; marg, comm, reg, endw, acts, endwm, endwms, endwf) = NamedTuple(Dict(Symbol(k) => sets[k] for k ∈ keys(sets)))

    fincome = NamedArray(mapslices(sum, hData["evfb"], dims=[1, 2])[1, 1, :], names(hData["evfb"])[[3]][1]) .- hData["vdep"]
    DPTAX = hData["vdpp"] .- hData["vdpb"]
    MPTAX = hData["vmpp"] .- hData["vmpb"]
    DGTAX = hData["vdgp"] .- hData["vdgb"]
    MGTAX = hData["vmgp"] .- hData["vmgb"]
    DITAX = hData["vdip"] .- hData["vdib"]
    MITAX = hData["vmip"] .- hData["vmib"]
    MFTAX = hData["vmfp"] .- hData["vmfb"]
    DFTAX = hData["vdfp"] .- hData["vdfb"]
    ETAX = hData["evfp"] .- hData["evfb"]
    PTAX = hData["makb"] .- hData["maks"]
    XTAXD = hData["vfob"] .- hData["vxsb"]
    MTAX = hData["vmsb"] .- hData["vcif"]

    TAXRIMP = mapslices(sum, MTAX, dims=[1, 2])[1, 1, :]
    TAXREXP = mapslices(sum, XTAXD, dims=[1, 3])[1, :, 1]
    TAXROUT = mapslices(sum, PTAX, dims=[1, 2])[1, 1, :]
    TAXRFU = mapslices(sum, ETAX, dims=[1, 2])[1, 1, :]
    TAXRIU = mapslices(sum, DFTAX .+ MFTAX, dims=[1, 2])[1, 1, :]
    TAXRIC = mapslices(sum, DITAX .+ MITAX, dims=1)[1, :]
    TAXRGC = mapslices(sum, DGTAX .+ MGTAX, dims=1)[1, :]
    TAXRPC = mapslices(sum, DPTAX .+ MPTAX, dims=1)[1, :]
    INDTAX = TAXRPC .+ TAXRGC .+ TAXRIC .+ TAXRIU .+ TAXRFU .+ TAXROUT .+ TAXREXP .+ TAXRIMP
    y = fincome .+ INDTAX

    yp = NamedArray(mapslices(sum, hData["vdpp"] .+ hData["vmpp"], dims=[1])[1, :], names(hData["vdpp"])[[2]][1])
    yg = NamedArray(mapslices(sum, hData["vdgp"] .+ hData["vmgp"], dims=[1])[1, :], names(hData["vdgp"])[[2]][1])
    walras_sup = sum(hData["vdib"]) + sum(hData["vmib"])
    walras_dem = walras_sup

    pint = NamedArray(ones(length(acts), length(reg)), (acts, reg))
    pva = NamedArray(ones(length(acts), length(reg)), (acts, reg))
    po = NamedArray(ones(length(acts), length(reg)), (acts, reg))
    pfa = NamedArray(ones(length(comm), length(acts), length(reg)), (comm, acts, reg))
    pfe = NamedArray(ones(length(endw), length(acts), length(reg)), (endw, acts, reg))
    pfd = NamedArray(ones(length(comm), length(acts), length(reg)), (comm, acts, reg))
    pfm = NamedArray(ones(length(comm), length(acts), length(reg)), (comm, acts, reg))
    ps = NamedArray(ones(length(comm), length(acts), length(reg)), (comm, acts, reg))
    pca = NamedArray(ones(length(comm), length(acts), length(reg)), (comm, acts, reg))
    pds = NamedArray(ones(length(comm), length(reg)), (comm, reg))
    peb = NamedArray(ones(length(endw), length(acts), length(reg)), (endw, acts, reg))
    ppm = NamedArray(ones(length(comm), length(reg)), (comm, reg))
    ppd = NamedArray(ones(length(comm), length(reg)), (comm, reg))
    ppa = NamedArray(ones(length(comm), length(reg)), (comm, reg))
    pga = NamedArray(ones(length(comm), length(reg)), (comm, reg))
    pes = NamedArray(ones(length(endw), length(acts), length(reg)), (endw, acts, reg))
    pgm = NamedArray(ones(length(comm), length(reg)), (comm, reg))
    pgd = NamedArray(ones(length(comm), length(reg)), (comm, reg))
    pga = NamedArray(ones(length(comm), length(reg)), (comm, reg))
    pid = NamedArray(ones(length(comm), length(reg)), (comm, reg))
    pim = NamedArray(ones(length(comm), length(reg)), (comm, reg))
    pia = NamedArray(ones(length(comm), length(reg)), (comm, reg))
    pms = NamedArray(ones(length(comm), length(reg)), (comm, reg))
    pcgdswld = 1
    psave = NamedArray(ones(length(reg)), (reg))

    up = NamedArray(100 * ones(length(reg)), (reg))

    pmds = NamedArray(ones(length(comm), length(reg), length(reg)), (comm, reg, reg))
    pcif = NamedArray(ones(length(comm), length(reg), length(reg)), (comm, reg, reg))
    pfob = NamedArray(ones(length(comm), length(reg), length(reg)), (comm, reg, reg))
    pgov = NamedArray(ones(length(reg)), (reg))
    pt = NamedArray(ones(length(marg)), (marg))
    pinv = NamedArray(ones(length(reg)), (reg))

    pe = NamedArray(ones(length(endwms), length(reg)), (endwms, reg))
    ptrans = NamedArray(ones(length(comm), length(reg), length(reg)), (comm, reg, reg))
    qe = NamedArray(mapslices(sum, hData["evos"], dims=[2])[:, 1, :], (endw, reg))[endwms, reg]

    pop = hData["pop"]
    ug = yg ./ pop

    #us = hData["save"] ./ pop

    #u = up .^ (yp ./ y) .* ug .^ (yg ./ y) .* us .^ (1 .- (yp .+ yg) ./ y)
    u = y ./ pop

    ppriv = NamedArray(ones(length(reg)), (reg))
    pfactor = NamedArray(ones(length(reg)), (reg))
    pfactwld = 1



    σ_vp = (hData["vdpp"] .+ hData["vmpp"]) ./ repeat(mapslices(sum, hData["vdpp"] .+ hData["vmpp"], dims=1), inner=[length(comm), 1])
    σ_vg = (hData["vdgp"] .+ hData["vmgp"]) ./ repeat(mapslices(sum, hData["vdgp"] .+ hData["vmgp"], dims=1), inner=[length(comm), 1])
    σ_vi = (hData["vdip"] .+ hData["vmip"]) ./ repeat(mapslices(sum, hData["vdip"] .+ hData["vmip"], dims=1), inner=[length(comm), 1])
    σ_vf = (hData["vdfp"] .+ hData["vmfp"]) ./ repeat(mapslices(sum, hData["vdfp"] .+ hData["vmfp"], dims=1), inner=[length(comm), 1, 1])
    σ_vdp = (hData["vdpp"]) ./ (hData["vdpp"] .+ hData["vmpp"])
    σ_vdg = (hData["vdgp"]) ./ (hData["vdgp"] .+ hData["vmgp"])
    σ_vdi = (hData["vdip"]) ./ (hData["vdip"] .+ hData["vmip"])
    σ_vdf = (hData["vdfp"]) ./ (hData["vdfp"] .+ hData["vmfp"])
    σ_qxs = (hData["vcif"]) ./ repeat(mapslices(sum, hData["vcif"], dims=2), inner=[1, length(reg), 1])
    σ_vff = (hData["evfp"]) ./ repeat(mapslices(sum, hData["evfp"], dims=1), inner=[length(endw), 1, 1])
    σsave = NamedArray(0.05 * ones(length(reg)), (reg))
    σ_qinv = NamedArray(0.05 * ones(length(reg)), (reg))

    for r ∈ reg
        σ_qinv[r] = hData["save"][r] / sum(hData["save"])
    end


    vfob = hData["vfob"]
    ϵ_qfa = NamedArray(ones(length(acts), length(reg)), (acts, reg))
    ϵ_qintva = NamedArray(ones(length(acts), length(reg)), (acts, reg))
    ϵ_qga = NamedArray(ones(length(reg)), (reg))
    ϵ_qia = NamedArray(ones(length(reg)), (reg))
    ϵ_qinv = 1
    σ_vtwr = NamedArray(0.01 * ones(length(marg), length(comm), length(reg), length(reg)), (marg, comm, reg, reg))
    for m ∈ marg
        for c ∈ comm
            for s ∈ reg
                for d ∈ reg
                    σ_vtwr[m, c, s, d] = hData["vtwr"][m, c, s, d] / hData["vcif"][c, s, d]
                end
            end
        end
    end

    σ_vif = NamedArray(0.5 * ones(length(acts), length(reg)), (acts, reg))
    σ_ρ = NamedArray(0.01 * ones(length(reg)), (reg))

    for r ∈ reg
        σ_ρ[r] = sum((hData["evos"][:, :, r])) / hData["vkb"][r]
    end

    rorg = 1
    rore = NamedArray(ones(length(reg)), reg)
    rorc = NamedArray(ones(length(reg)), reg)
    rental = NamedArray(ones(length(reg)), reg)

    p = NamedArray(ones(length(reg)), reg)

    return (data = Dict(
        "ppriv" => ppriv,
        "pfactor" => pfactor,
        "pint" => pint,
        "pva" => pva,
        "po" => po,
        "pfa" => pfa,
        "pfe" => pfe,
        "pfm" => pfm,
        "pfd" => pfd,
        "ps" => ps,
        "pca" => pca,
        "pds" => pds,
        "peb" => peb,
        "ppm" => ppm,
        "ppd" => ppd,
        "ppa" => ppa,
        "pga" => pga,
        "pes" => pes,
        "pgm" => pgm,
        "pgd" => pgd,
        "pga" => pga,
        "pid" => pid,
        "pim" => pim,
        "pia" => pia,
        "pms" => pms,
        "pcgdswld" => pcgdswld,
        "psave" => psave,
        "u" => u,
        "up" => up,
        "ug" => ug,
        #"us" => us,
        "pmds" => pmds,
        "pcif" => pcif,
        "pfob" => pfob,
        "pgov" => pgov,
        "pt" => pt,
        "pinv" => pinv,
        "pe" => pe,
        "ptrans" => ptrans,
        "fincome" => fincome,
        "y" => y,
        "yp" => yp,
        "yg" => yg,
        "walras_sup" => walras_sup,
        "walras_dem" => walras_dem,
        "qe" => qe,
        "pop" => pop,
        "pfactwld" => pfactwld,
        "σ_vp" => σ_vp,
        "σ_vg" => σ_vg,
        "σ_vi" => σ_vi,
        "σ_vf" => σ_vf,
        "σ_vdp" => σ_vdp,
        "σ_vdg" => σ_vdg,
        "σ_vdi" => σ_vdi,
        "σ_vdf" => σ_vdf,
        "σ_qxs" => σ_qxs,
        "σ_vff" => σ_vff,
        "σsave" => σsave,
        "σ_qinv" => σ_qinv,
        "vfob" => vfob,
        "ϵ_qfa" => ϵ_qfa,
        "ϵ_qintva" => ϵ_qintva,
        "ϵ_qga" => ϵ_qga,
        "ϵ_qia" => ϵ_qia,
        "ϵ_qinv" => ϵ_qinv,
        "σ_vtwr" => σ_vtwr,
        "σ_vif" => σ_vif,
        "σ_ρ" => σ_ρ,
        "rorg" => rorg,
        "rore" => rore,
        "rorc" => rorc,
        "rental" => rental,
        "p" => p
    ))

end