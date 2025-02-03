# This calculates some initial values 
function prepare_initial_values(; sets, hData, hParameters)
    # Read the sets
    (; marg, comm, reg, endw, acts, endwm, endwms, endwf) = NamedTuple(Dict(Symbol(k) => sets[k] for k ∈ keys(sets)))

    fincome = NamedArray(mapslices(sum, hData["evfb"], dims=[1, 2])[1, 1, :], names(hData["evfb"])[[3]][1]) .- hData["vdep"]
    DPTAX = hData["vdpp"]  .- hData["vdpb"]
    MPTAX = hData["vmpp"]  .- hData["vmpb"]
    DGTAX = hData["vdgp"]  .- hData["vdgb"]
    MGTAX = hData["vmgp"]  .- hData["vmgb"]
    DITAX = hData["vdip"]  .- hData["vdib"]
    MITAX = hData["vmip"]  .- hData["vmib"]
    MFTAX = hData["vmfp"] .- hData["vmfb"]
    DFTAX = hData["vdfp"] .- hData["vdfb"]
    ETAX = hData["evfp"] .- hData["evfb"]
    PTAX = hData["makb"] .- hData["maks"]
    XTAXD = hData["vfob"] .- hData["vxsb"]
    MTAX = hData["vmsb"] .- hData["vcif"]
    
    TAXRIMP = mapslices(sum, MTAX, dims = [1,2])[1,1,:]
    TAXREXP = mapslices(sum, XTAXD, dims = [1,3])[1,:,1]
    TAXROUT = mapslices(sum, PTAX, dims = [1,2])[1,1,:]
    TAXRFU = mapslices(sum, ETAX, dims = [1,2])[1,1,:]
    TAXRIU = mapslices(sum, DFTAX .+ MFTAX, dims = [1,2])[1,1,:]
    TAXRIC = mapslices(sum, DITAX .+ MITAX, dims = 1)[1,:]
    TAXRGC = mapslices(sum, DGTAX .+ MGTAX, dims = 1)[1,:]
    TAXRPC = mapslices(sum, DPTAX .+ MPTAX, dims = 1)[1,:]
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
    u = NamedArray(ones(length(reg)), (reg))
    up = NamedArray(ones(length(reg)), (reg))
    ug = NamedArray(ones(length(reg)), (reg))
    us = NamedArray(ones(length(reg)), (reg))
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

    ppriv = NamedArray(ones(length(reg)), (reg))
    pfactor = NamedArray(ones(length(reg)), (reg))
    pfactwld = 1


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
        "us" => us,
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
        "pfactwld" => pfactwld
    ))

end