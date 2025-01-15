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
        end for (k) in ["to", "tfe", "tx", "txs", "tm", "tms", "tfd", "tfm", "tpd", "tpm", "tgd", "tgm", "tid", "tim", "tinc", "qesf", "qe", "ppa", "α_qintva", "γ_qintva", "α_qfa", "γ_qfa", "α_qfe", "γ_qfe", "α_qfdqfm", "γ_qfdqfm", "α_qca", "γ_qca", "α_pca", "γ_pca", "σyp", "σyg", "β_qpa", "α_qpdqpm", "γ_qpdqpm", "α_qga", "γ_qga", "α_qgdqgm", "γ_qgdqgm", "α_qia", "γ_qia", "α_qidqim", "γ_qidqim", "α_qxs", "γ_qxs", "α_qtmfsd", "α_qst", "γ_qst", "α_qes2", "γ_qes2", "α_qinv", "δ", "ρ", "pop"]
    )
    
    fixed["ppa"][:, :] .= false
    ## The price of the first commodity in the first region is fixed
    fixed["ppa"][comm[1], reg[1]] = true

    # Calculate calibrated_parameter initial values
    return (sets=sets, parameters=parameters, data=data, fixed=fixed)

end