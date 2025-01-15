function aggregate_data(; hData, hParameters, hSets, comMap, regMap, endMap)

        marMap = comMap[map(f -> in(f, hSets["marg"]), names(comMap)[1])]

        # Fixed map
        fixedMap = NamedArray(["mobile", "sluggish", "fixed"],["mobile", "sluggish", "fixed"])

        # Generate the aggregated hData
        dataAg =
                Dict(
                        "vdpb" => GeneralEquilibrium.agg(hData["vdpb"], [comMap, regMap]),
                        "vmpb" => GeneralEquilibrium.agg(hData["vmpb"], [comMap, regMap]),
                        "vdpp" => GeneralEquilibrium.agg(hData["vdpp"], [comMap, regMap]),
                        "vmpp" => GeneralEquilibrium.agg(hData["vmpp"], [comMap, regMap]),
                        "vdgb" => GeneralEquilibrium.agg(hData["vdgb"], [comMap, regMap]),
                        "vmgb" => GeneralEquilibrium.agg(hData["vmgb"], [comMap, regMap]),
                        "vdgp" => GeneralEquilibrium.agg(hData["vdgp"], [comMap, regMap]),
                        "vmgp" => GeneralEquilibrium.agg(hData["vmgp"], [comMap, regMap]),
                        "vdib" => GeneralEquilibrium.agg(hData["vdib"], [comMap, regMap]),
                        "vmib" => GeneralEquilibrium.agg(hData["vmib"], [comMap, regMap]),
                        "vdip" => GeneralEquilibrium.agg(hData["vdip"], [comMap, regMap]),
                        "vmip" => GeneralEquilibrium.agg(hData["vmip"], [comMap, regMap]),
                        "vdfb" => GeneralEquilibrium.agg(hData["vdfb"], [comMap, comMap, regMap]),
                        "vmfb" => GeneralEquilibrium.agg(hData["vmfb"], [comMap, comMap, regMap]),
                        "vdfp" => GeneralEquilibrium.agg(hData["vdfp"], [comMap, comMap, regMap]),
                        "vmfp" => GeneralEquilibrium.agg(hData["vmfp"], [comMap, comMap, regMap]),
                        "evfb" => GeneralEquilibrium.agg(hData["evfb"], [endMap, comMap, regMap]),
                        "evfp" => GeneralEquilibrium.agg(hData["evfp"], [endMap, comMap, regMap]),
                        "evos" => GeneralEquilibrium.agg(hData["evos"], [endMap, comMap, regMap]),
                        "vmsb" => GeneralEquilibrium.agg(hData["vmsb"], [comMap, regMap, regMap]),
                        "vxsb" => GeneralEquilibrium.agg(hData["vxsb"], [comMap, regMap, regMap]),
                        "vfob" => GeneralEquilibrium.agg(hData["vfob"], [comMap, regMap, regMap]),
                        "vcif" => GeneralEquilibrium.agg(hData["vcif"], [comMap, regMap, regMap]),
                        "vtwr" => GeneralEquilibrium.agg(hData["vtwr"], [marMap, comMap, regMap, regMap]),
                        "vst" => GeneralEquilibrium.agg(hData["vst"], [marMap, regMap]),
                        "save" => GeneralEquilibrium.agg(hData["save"], [regMap]),
                        "vdep" => GeneralEquilibrium.agg(hData["vdep"], [regMap]),
                        "vkb" => GeneralEquilibrium.agg(hData["vkb"], [regMap]),
                        "maks" => GeneralEquilibrium.agg(hData["maks"], [comMap, comMap, regMap]),
                        "makb" => GeneralEquilibrium.agg(hData["makb"], [comMap, comMap, regMap]),
                        "pop" => GeneralEquilibrium.agg(hData["pop"], [regMap])
                        )

        # Generate the aggregated hParameters
        paramAg = Dict(
                ### "esbv" => GeneralEquilibrium.aggComb(hParameters["esbv"], NamedArray(mapslices(sum, hData["evfp"], dims=1)[1, :, :], names(hParameters["esbv"])), [comMap, regMap]),
                "esbv" => GeneralEquilibrium.aggComb(hParameters["esbv"], NamedArray(repeat(reshape(mapslices(sum, hData["evfp"], dims=[1, 3])[1, :, 1], (length(comMap), 1)), inner=[1, length(regMap)]), names(hParameters["esbv"])), [comMap, regMap]),
                "esbt" => GeneralEquilibrium.aggComb(hParameters["esbt"], NamedArray(
                                mapslices(sum, hData["vdfp"], dims=1)[1, :, :]
                                + mapslices(sum, hData["vmfp"], dims=1)[1, :, :]
                                + mapslices(sum, hData["evfp"], dims=1)[1, :, :], names(hParameters["esbt"])), [comMap, regMap]),
                "esbc" => GeneralEquilibrium.aggComb(hParameters["esbt"], NamedArray(
                                mapslices(sum, hData["vdfp"], dims=1)[1, :, :] + mapslices(sum, hData["vmfp"], dims=1)[1, :, :], names(hParameters["esbv"])), [comMap, regMap]),
                "etrq" => GeneralEquilibrium.aggComb(hParameters["etrq"], NamedArray(
                                mapslices(sum, hData["maks"], dims=1)[1, :, :], names(hParameters["etrq"])), [comMap, regMap]),
                "esbq" => GeneralEquilibrium.aggComb(hParameters["esbq"], NamedArray(
                                mapslices(sum, hData["maks"], dims=2)[:, 1, :], names(hParameters["esbq"])), [comMap, regMap]),
                "esbg" => GeneralEquilibrium.aggComb(hParameters["esbg"], NamedArray(mapslices(sum, hData["vdgp"] .+ hData["vmgp"], dims=1)[1, :], names(hParameters["esbg"])[1]), [regMap]),
                #"esbd" => GeneralEquilibrium.aggComb(hParameters["esbd"], NamedArray(mapslices(sum, hData["vdfp"] .+ hData["vmfp"], dims=2)[:, 1, :], names(hParameters["esbd"])) .+ hData["vdpp"].+ hData["vmpp"] .+ hData["vdgp"].+ hData["vmgp"] .+ hData["vdip"].+ hData["vmip"], [comMap, regMap]),
                #"esbm" => GeneralEquilibrium.aggComb(hParameters["esbm"], NamedArray(mapslices(sum, hData["vmfp"], dims=2)[:, 1, :], names(hParameters["esbm"])) .+ hData["vmpp"] .+ hData["vmgp"] .+ hData["vmip"], [comMap, regMap]),
                "esbd" => GeneralEquilibrium.aggComb(hParameters["esbd"], NamedArray(repeat(reshape(mapslices(sum, NamedArray(mapslices(sum, hData["vdfp"] .+ hData["vmfp"], dims=2)[:, 1, :], names(hParameters["esbd"])) .+ hData["vdpp"] .+ hData["vmpp"] .+ hData["vdgp"] .+ hData["vmgp"] .+ hData["vdip"] .+ hData["vmip"], dims=2), (length(comMap), 1)), inner=[1, length(regMap)]), names(hParameters["esbd"])), [comMap, regMap]),
                "esbm" => GeneralEquilibrium.aggComb(hParameters["esbm"], NamedArray(repeat(reshape(mapslices(sum, NamedArray(mapslices(sum, hData["vmfp"], dims=2)[:, 1, :], names(hParameters["esbm"])) .+ hData["vmpp"] .+ hData["vmgp"] .+ hData["vmip"], dims=2), (length(comMap), 1)), inner=[1, length(regMap)]), names(hParameters["esbm"])), [comMap, regMap]),
                "esbs" => GeneralEquilibrium.aggComb(hParameters["esbs"], NamedArray(mapslices(sum, hData["vst"], dims=2)[:, 1], names(hParameters["esbs"])[1]), [marMap]),
                "subp" => GeneralEquilibrium.aggComb(hParameters["subp"], hData["vdpp"] .+ hData["vmpp"], [comMap, regMap]),
                "incp" => GeneralEquilibrium.aggComb(hParameters["incp"], hData["vdpp"] .+ hData["vmpp"], [comMap, regMap]),
                "rflx" => hParameters["rflx"],
                "etre" => GeneralEquilibrium.aggComb(hParameters["etre"], NamedArray(mapslices(sum, hData["evos"], dims=2)[:, 1, :], names(hParameters["etre"])), [endMap, regMap]), #hParameters["etre"],
                "eflg" => GeneralEquilibrium.agg(hParameters["eflg"], [endMap, fixedMap])
        )

        paramAg["eflg"] = paramAg["eflg"] ./ maximum.([paramAg["eflg"][i, :] for i ∈ 1:size(paramAg["eflg"], 1)])
        paramAg["eflg"][paramAg["eflg"].<0] .= 0


        reg = unique(regMap)
        comm = unique(comMap)
        marg = unique(marMap)
        acts = unique(comMap)
        endw = unique(endMap)


        # Return the aggregated hData
        return (hData=dataAg, hParameters=paramAg, hSets=Dict("reg" => reg, "comm" => comm, "marg" => marg, "acts" => acts, "endw" => endw))
end
