function generate_starting_values(; hSets, hData, hParameters)

    (sets) = prepare_sets(hSets=hSets, hParameters=hParameters)

    (parameters) = prepare_parameters(hParameters=hParameters)

    data = prepare_initial_values(sets=sets, hData=hData, hParameters=hParameters)

    data = prepare_taxes(data=data, hData=hData)

    data = prepare_quantities(data=data, parameters=parameters, sets=sets, hData=hData)

    (; parameters, data) = prepare_initial_calibrated_parameters(data=data, sets=sets, parameters=parameters, hData=hData)

    # Prepare a set of fixed parameters
    (; comm, reg) = NamedTuple(Dict(Symbol(k) => sets[k] for k ∈ ["comm", "reg"]))
    fixed = Dict(
        k => begin
            if data[k] isa AbstractVector
                NamedArray(trues(size(data[k])), names(data[k])...)
            elseif data[k] isa AbstractArray
                NamedArray(trues(size(data[k])), names(data[k]))
            else
                true
            end
        end for (k) in ["to", "tfe", "tx", "txs", "tm", "tms", "tfd", "tfm", "tpd", "tpm", "tgd", "tgm", "tid", "tim", "tinc", "qesf", "qe", "α_qintva", "γ_qintva", "α_qfa", "γ_qfa", "α_qfe", "γ_qfe", "α_qfdqfm", "γ_qfdqfm", "α_qca", "γ_qca", "α_pca", "γ_pca", "σyp", "σyg", "β_qpa", "α_qpdqpm", "γ_qpdqpm", "α_qga", "γ_qga", "α_qgdqgm", "γ_qgdqgm", "α_qia", "γ_qia", "α_qidqim", "γ_qidqim", "α_qxs", "γ_qxs", "α_qtmfsd", "α_qst", "γ_qst", "α_qes2", "γ_qes2", "α_qinv", "δ", "ρ", "pop"]
    )

    fixed = merge(Dict("pfactwld"=>true), fixed)

    #fixed["ppa"][:, :] .= false
    ## The price of the first commodity in the first region is fixed
    #fixed["ppa"][comm[1], reg[1]] = true

    q_min = 1e-8
    q_max = 1e+12
    p_min = 1e-8
    p_max = 1e+8
    t_min = 1e-8
    t_max = 1e+2
    α_min = 1e-8
    α_max = 1
    γ_min = 1e-8
    ϵ_min = 1e-8
    σ_min = 0
    y_min = 1e-8
    y_max = 1e12
    v_min = 0


    lower = merge(
        Dict(k => q_min for k ∈ ["pop", "qint", "qva", "qo", "qfa", "qfe", "qfd", "qfm", "qca", "qc", "qes", "u", "up", "ug", "us", "qpa", "qpd", "qpm", "qga", "qgd", "qgm", "qia", "qid", "qim", "qinv", "qms", "qxs", "qtmfsd", "qtm", "qst", "qds", "qe", "qesf", "walras_sup", "walras_dem", "kb", "ke", "globalcgds"]),
        Dict("qsave" => -q_max),
        Dict(k => p_min for k ∈ ["pfactwld","pca", "pcgdswld", "pcif", "pds", "pe", "peb", "pes", "pfa", "pfactor", "pfd", "pfe", "pfm", "pfob", "pga", "pgd", "pgm", "pgov", "pia", "pid", "pim", "pint", "pinv", "pmds", "pms", "po", "pop", "ppa", "ppd", "ppm", "ppriv", "ps", "psave", "pt", "ptrans", "pva"]),
        Dict(k => t_min for k ∈ ["tfd", "tfe", "tfm", "tgd", "tgm", "tid", "tim", "tinc", "tm", "tms", "to", "tpd", "tpm", "tx", "txs"]),
        Dict(k => α_min for k ∈ ["α_pca", "α_qca", "α_qes2", "α_qfa", "α_qfdqfm", "α_qfe", "α_qga", "α_qgdqgm", "α_qia", "α_qidqim", "α_qintva", "α_qinv", "α_qpdqpm", "α_qst", "α_qtmfsd", "α_qxs"]),
        Dict(k => γ_min for k ∈ ["γ_pca", "γ_qca", "γ_qes2", "γ_qfa", "γ_qfdqfm", "γ_qfe", "γ_qga", "γ_qgdqgm", "γ_qia", "γ_qidqim", "γ_qintva", "γ_qpdqpm", "γ_qst", "γ_qxs"]),
        Dict(k => ϵ_min for k ∈ ["ϵ_qes2", "ϵ_qfa", "ϵ_qfdqfm", "ϵ_qfe", "ϵ_qga", "ϵ_qgdqgm", "ϵ_qia", "ϵ_qidqim", "ϵ_qintva", "ϵ_qinv", "ϵ_qpdqpm", "ϵ_qxs"]),
        Dict(k => σ_min for k ∈ ["σ_qxs", "σ_vdf", "σ_vdg", "σ_vdi", "σ_vdp", "σ_vf", "σ_vff", "σ_vg", "σ_vi", "σ_vif", "σ_vp", "σ_vtwr", "σ_ρ", "σsave", "σyg", "σyp"]),
        Dict("σ_qinv" => -1),
        Dict(k => y_min for k ∈ ["fincome", "y", "yp", "yg"]),
        Dict(k => v_min for k ∈ ["vdfp", "vmfp", "vdpp", "vmpp", "vdgp", "vmgp", "vdip", "vmip", "evfp", "evos", "vfob", "vcif", "vst", "vtwr", "maks", "vkb"]),
        Dict(k => 0 for k ∈ ["uelas", "uepriv"]),
        Dict("β_qpa" => 1e-8),
        Dict(k => 0 for k ∈ ["δ", "ρ"])
    )

    upper = merge(
        Dict(k => q_max for k ∈ ["pop", "qint", "qva", "qo", "qfa", "qfe", "qfd", "qfm", "qca", "qc", "qes", "u", "up", "ug", "us", "qpa", "qpd", "qpm", "qga", "qgd", "qgm", "qsave", "qia", "qid", "qim", "qinv", "qms", "qxs", "qtmfsd", "qtm", "qst", "qds", "qe", "qesf", "walras_sup", "walras_dem", "kb", "ke", "globalcgds"]),
        Dict(k => p_max for k ∈ ["pfactwld","pca", "pcgdswld", "pcif", "pds", "pe", "peb", "pes", "pfa", "pfactor", "pfd", "pfe", "pfm", "pfob", "pga", "pgd", "pgm", "pgov", "pia", "pid", "pim", "pint", "pinv", "pmds", "pms", "po", "pop", "ppa", "ppd", "ppm", "ppriv", "ps", "psave", "pt", "ptrans", "pva"]),
        Dict(k => t_max for k ∈ ["tfd", "tfe", "tfm", "tgd", "tgm", "tid", "tim", "tinc", "tm", "tms", "to", "tpd", "tpm", "tx", "txs"]),
        Dict(k => α_max for k ∈ ["α_pca", "α_qca", "α_qes2", "α_qfa", "α_qfdqfm", "α_qfe", "α_qga", "α_qgdqgm", "α_qia", "α_qidqim", "α_qintva", "α_qinv", "α_qpdqpm", "α_qst", "α_qtmfsd", "α_qxs"]),
        Dict("σ_qinv" => 1),
        Dict(k => y_max for k ∈ ["fincome", "y", "yp", "yg"]),
        Dict(k => 10 for k ∈ ["uelas", "uepriv"]),
        Dict(k => 1 for k ∈ ["δ", "ρ"])
    )


    # Calculate calibrated_parameter initial values
    return (sets=sets, parameters=parameters, data=data, fixed=fixed, lower=lower, upper=upper)

end