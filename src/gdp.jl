function gdp(; sets, data0, data1, parameters, max_iter=50, constr_viol_tol=1e-5, bound_push=1e-15)


    return NamedArray((mapslices(sum, data1["qga"] .* data0["pga"] .+ data1["qpa"] .* data0["ppa"] .+ data1["qia"] .* data0["pia"] .+
                mapslices(sum, data1["qxs"] .* data0["pfob"], dims=2)[:, 1, :] .+
                mapslices(sum, data1["qxs"] .* data0["pcif"], dims=3)[:, :, 1], dims=1)[1, :] .+ data1["pds"]["svces", :] .* data0["qst"][1, :]
) ./(mapslices(sum, data0["qga"] .* data0["pga"] .+ data0["qpa"] .* data0["ppa"] .+ data0["qia"] .* data0["pia"] .+
mapslices(sum, data0["qxs"] .* data0["pfob"], dims=2)[:, 1, :] .+
mapslices(sum, data0["qxs"] .* data0["pcif"], dims=3)[:, :, 1], dims=1)[1, :] .+ data0["pds"]["svces", :] .* data0["qst"][1, :]
),sets["reg"])

    
end
