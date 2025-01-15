function prepare_parameters(; hParameters)
    # Prepare parameters
    esubt = hParameters["esbt"]
    esubc = hParameters["esbc"]
    esubva = hParameters["esbv"]
    esubd = hParameters["esbd"]
    etraq = hParameters["etrq"]
    esubq = hParameters["esbq"]
    subpar = hParameters["subp"]
    incpar = hParameters["incp"]
    esubg = hParameters["esbg"]
    esubm = hParameters["esbm"]
    esubs = hParameters["esbs"]
    etrae = hParameters["etre"]
    endowflag = hParameters["eflg"]

    return (parameters = Dict(
        "esubt" => esubt,
        "esubc" => esubc,
        "esubva" => esubva,
        "esubd" => esubd,
        "etraq" => etraq,
        "esubq" => esubq,
        "subpar" => subpar,
        "incpar" => incpar,
        "esubg" => esubg,
        "esubm" => esubm,
        "esubs" => esubs,
        "etrae" => etrae,
        "endowflag" => endowflag
    ))

end