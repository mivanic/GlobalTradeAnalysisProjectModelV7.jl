"""
        aggregate_data(; hData, hParameters, hSets, comMap, regMap, endMap)

Aggregates data, parameters and sets based on the provided mapping vectors.

        Args:
                hData: a dictionary of GTAP data (arrays) with the names as found in HAR files
                hParameters: a dictionary of GTAP parameters (arrays) with the names as found in HAR files
                hSets: a dictionary of GTAP sets (vectors of strings) with the names as found in HAR files
                comMap: a mapping vector (NamedVector) which provides for each element in set `comm` the desired aggregate name
                regMap: a mapping vector (NamedVector) which provides for each element in set `reg` the desired aggregate name
                endMap: a mapping vector (NamedVector) which provides for each element in set `endw` the desired aggregate name

        Retruns:
                A named tuple with elements hData, hParameters, hSets containing the aggregated data
        
        Example:

        ```julia
        (; hData, hParameters, hSets) = aggregate_data(hData=data, hParameters=parameters, hSets=sets, comMap=comMap, regMap=regMap, endMap=endMap)
        ```

"""
function aggregate_data(; hData, hParameters, hSets, comMap, regMap, endMap)

        marMap = comMap[map(f -> in(f, hSets["marg"]), names(comMap)[1])]

        # Fixed map
        fixedMap = NamedArray(["mobile", "sluggish", "fixed"],["mobile", "sluggish", "fixed"])

        # Generate the aggregated hData
        dataAg =
                Dict(
                        "vdpb" => agg(hData["vdpb"], [comMap, regMap]),
                        "vmpb" => agg(hData["vmpb"], [comMap, regMap]),
                        "vdpp" => agg(hData["vdpp"], [comMap, regMap]),
                        "vmpp" => agg(hData["vmpp"], [comMap, regMap]),
                        "vdgb" => agg(hData["vdgb"], [comMap, regMap]),
                        "vmgb" => agg(hData["vmgb"], [comMap, regMap]),
                        "vdgp" => agg(hData["vdgp"], [comMap, regMap]),
                        "vmgp" => agg(hData["vmgp"], [comMap, regMap]),
                        "vdib" => agg(hData["vdib"], [comMap, regMap]),
                        "vmib" => agg(hData["vmib"], [comMap, regMap]),
                        "vdip" => agg(hData["vdip"], [comMap, regMap]),
                        "vmip" => agg(hData["vmip"], [comMap, regMap]),
                        "vdfb" => agg(hData["vdfb"], [comMap, comMap, regMap]),
                        "vmfb" => agg(hData["vmfb"], [comMap, comMap, regMap]),
                        "vdfp" => agg(hData["vdfp"], [comMap, comMap, regMap]),
                        "vmfp" => agg(hData["vmfp"], [comMap, comMap, regMap]),
                        "evfb" => agg(hData["evfb"], [endMap, comMap, regMap]),
                        "evfp" => agg(hData["evfp"], [endMap, comMap, regMap]),
                        "evos" => agg(hData["evos"], [endMap, comMap, regMap]),
                        "vmsb" => agg(hData["vmsb"], [comMap, regMap, regMap]),
                        "vxsb" => agg(hData["vxsb"], [comMap, regMap, regMap]),
                        "vfob" => agg(hData["vfob"], [comMap, regMap, regMap]),
                        "vcif" => agg(hData["vcif"], [comMap, regMap, regMap]),
                        "vtwr" => agg(hData["vtwr"], [marMap, comMap, regMap, regMap]),
                        "vst" => agg(hData["vst"], [marMap, regMap]),
                        "save" => agg(hData["save"], [regMap]),
                        "vdep" => agg(hData["vdep"], [regMap]),
                        "vkb" => agg(hData["vkb"], [regMap]),
                        "maks" => agg(hData["maks"], [comMap, comMap, regMap]),
                        "makb" => agg(hData["makb"], [comMap, comMap, regMap]),
                        "pop" => agg(hData["pop"], [regMap])
                        )

        # Generate the aggregated hParameters
        paramAg = Dict(
                "esbv" => aggComb(hParameters["esbv"], NamedArray(mapslices(sum, hData["evfp"], dims=1)[1, :, :], names(hParameters["esbv"])), [comMap, regMap]),
                #"esbv" => aggComb(hParameters["esbv"], NamedArray(repeat(reshape(mapslices(sum, hData["evfp"], dims=[1, 3])[1, :, 1], (length(comMap), 1)), inner=[1, length(regMap)]), names(hParameters["esbv"])), [comMap, regMap]),
                "esbt" => aggComb(hParameters["esbt"], NamedArray(
                                mapslices(sum, hData["vdfp"], dims=1)[1, :, :]
                                + mapslices(sum, hData["vmfp"], dims=1)[1, :, :]
                                + mapslices(sum, hData["evfp"], dims=1)[1, :, :], names(hParameters["esbt"])), [comMap, regMap]),
                "esbc" => aggComb(hParameters["esbt"], NamedArray(
                                mapslices(sum, hData["vdfp"], dims=1)[1, :, :] + mapslices(sum, hData["vmfp"], dims=1)[1, :, :], names(hParameters["esbv"])), [comMap, regMap]),
                "etrq" => aggComb(hParameters["etrq"], NamedArray(
                                mapslices(sum, hData["maks"], dims=1)[1, :, :], names(hParameters["etrq"])), [comMap, regMap]),
                "esbq" => aggComb(hParameters["esbq"], NamedArray(
                                mapslices(sum, hData["maks"], dims=2)[:, 1, :], names(hParameters["esbq"])), [comMap, regMap]),
                "esbg" => aggComb(hParameters["esbg"], NamedArray(mapslices(sum, hData["vdgp"] .+ hData["vmgp"], dims=1)[1, :], names(hParameters["esbg"])[1]), [regMap]),
                "esbd" => aggComb(hParameters["esbd"], NamedArray(mapslices(sum, hData["vdfp"] .+ hData["vmfp"], dims=2)[:, 1, :], names(hParameters["esbd"])) .+ hData["vdpp"].+ hData["vmpp"] .+ hData["vdgp"].+ hData["vmgp"] .+ hData["vdip"].+ hData["vmip"], [comMap, regMap]),
                "esbm" => aggComb(hParameters["esbm"], NamedArray(mapslices(sum, hData["vmfp"], dims=2)[:, 1, :], names(hParameters["esbm"])) .+ hData["vmpp"] .+ hData["vmgp"] .+ hData["vmip"], [comMap, regMap]),
                #"esbd" => aggComb(hParameters["esbd"], NamedArray(repeat(reshape(mapslices(sum, NamedArray(mapslices(sum, hData["vdfp"] .+ hData["vmfp"], dims=2)[:, 1, :], names(hParameters["esbd"])) .+ hData["vdpp"] .+ hData["vmpp"] .+ hData["vdgp"] .+ hData["vmgp"] .+ hData["vdip"] .+ hData["vmip"], dims=2), (length(comMap), 1)), inner=[1, length(regMap)]), names(hParameters["esbd"])), [comMap, regMap]),
                #"esbm" => aggComb(hParameters["esbm"], NamedArray(repeat(reshape(mapslices(sum, NamedArray(mapslices(sum, hData["vmfp"], dims=2)[:, 1, :], names(hParameters["esbm"])) .+ hData["vmpp"] .+ hData["vmgp"] .+ hData["vmip"], dims=2), (length(comMap), 1)), inner=[1, length(regMap)]), names(hParameters["esbm"])), [comMap, regMap]),
                "esbs" => aggComb(hParameters["esbs"], NamedArray(mapslices(sum, hData["vst"], dims=2)[:, 1], names(hParameters["esbs"])[1]), [marMap]),
                "subp" => aggComb(hParameters["subp"], hData["vdpp"] .+ hData["vmpp"], [comMap, regMap]),
                "incp" => aggComb(hParameters["incp"], hData["vdpp"] .+ hData["vmpp"], [comMap, regMap]),
                "rflx" => hParameters["rflx"],
                "etre" => aggComb(hParameters["etre"], NamedArray(mapslices(sum, hData["evos"], dims=2)[:, 1, :], names(hParameters["etre"])), [endMap, regMap]), #hParameters["etre"],
                "eflg" => agg(hParameters["eflg"], [endMap, fixedMap])
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
