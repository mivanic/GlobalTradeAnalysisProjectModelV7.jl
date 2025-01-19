function calculate_gdp(; sets, data0, data1)

    return NamedArray((mapslices(sum, data1["qga"] .* data0["pga"] .+ data1["qpa"] .* data0["ppa"] .+ data1["qia"] .* data0["pia"] .+
                                      mapslices(sum, ifelse.(is.nan.(data1["qxs"] .* data0["pfob"]), 0.0, data1["qxs"] .* data0["pfob"]), dims=2)[:, 1, :] .+
                                      mapslices(sum, ifelse.(is.nan.(data1["qxs"] .* data0["pcif"]), 0.0, data1["qxs"] .* data0["pcif"]), dims=3)[:, :, 1], dims=1)[1, :] .+ data1["pds"]["svces", :] .* data0["qst"][1, :]
        ), sets["reg"])


end
