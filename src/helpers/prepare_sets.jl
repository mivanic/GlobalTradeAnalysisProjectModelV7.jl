function prepare_sets(; hSets, hParameters)
    # generate all sets
    (; marg, comm, reg, endw, acts) = NamedTuple(Dict(Symbol(k) => hSets[k] for k ∈ keys(hSets)))

    endwc = ["capital"]
    endws = names(hParameters["eflg"][hParameters["eflg"][:, "sluggish"].==1, :])[1]
    endwm = names(hParameters["eflg"][hParameters["eflg"][:, "mobile"].==1, :])[1]
    endwms = [endwm; endws]
    endwf = names(hParameters["eflg"][hParameters["eflg"][:, "fixed"].==1, :])[1]

    return (sets = Dict("endwc" => endwc, "endws" => endws, "endwm" => endwm, "endwms" => endwms, "endwf" => endwf, "marg" => marg, "comm" => comm, "reg" => reg, "acts" => acts, "endw" => endw))
end