
"""
    calculate_gdp(; sets, data0, data1)
"""
function calculate_gdp(; sets, data0, data1)

    return NamedArray(
        (
            mapslices(sum,
                data1["qga"] .* data0["pga"] .+ data1["qpa"] .* data0["ppa"] .+ data1["qia"] .* data0["pia"] .+
                mapslices(sum, ifelse.(isnan.(data1["qxs"] .* data0["pfob"]), 0.0, data1["qxs"] .* data0["pfob"]), dims=3)[:, :, 1] .-
                mapslices(sum, ifelse.(isnan.(data1["qxs"] .* data0["pcif"]), 0.0, data1["qxs"] .* data0["pcif"]), dims=2)[:, 1, :], dims=1)[1, :] .+ mapslices(sum, data0["pds"][sets["marg"], :] .* data1["qst"][sets["marg"], :], dims = 1)[1,:]
        ), sets["reg"])
end
