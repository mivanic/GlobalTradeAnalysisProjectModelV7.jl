function generate_calibration_inputs(model_container, start_data)

    # We do not want to mess with user's fixed dictionary
    fixed = deepcopy(model_container.fixed)
    calibrate_start = deepcopy(model_container.data)
    data = deepcopy(model_container.data)
    parameters = deepcopy(model_container.parameters)
    sets = deepcopy(model_container.sets)

    # CAL-I
    fixed["α_qxs"] .= false
    fixed["σ_qxs"] = NamedArray(trues(size(data["σ_qxs"])), names(data["σ_qxs"]))
    for r ∈ sets["reg"]
        for c ∈ sets["comm"]
            mr = findall(.!isnan.(data["vcif"][c, :, r]))[1]
            fixed["σ_qxs"][c, mr, r] .= false
        end
    end
    fixed["ϵ_qxs"] = NamedArray(trues(size(data["ϵ_qxs"])), names(data["ϵ_qxs"]))

    calibrate_start["σ_qxs"] = (start_data["vcif"]) ./ repeat(reshape(mapslices(sum, (start_data["vcif"]), dims=2)[:, 1, :], [size(start_data["vcif"])[1], 1, size(start_data["vcif"])[3]]...), inner=[1, size(start_data["vcif"])[2], 1])

    # CAL-II
    fixed["α_qintva"] .= false
    fixed["α_qfe"] .= false
    fixed["σ_vff"] = NamedArray(trues(size(data["σ_vff"])), names(data["σ_vff"]))
    fixed["σ_vff"]["capital", :, :] .= false
    fixed["σ_vif"] = NamedArray(trues(size(data["σ_vif"])), names(data["σ_vif"]))
    fixed["ϵ_qintva"] = NamedArray(trues(size(data["ϵ_qintva"])), names(data["ϵ_qintva"]))
    fixed["ϵ_qfe"] = NamedArray(trues(size(data["ϵ_qfe"])), names(data["ϵ_qfe"]))

    calibrate_start["σ_vff"] = (start_data["evfp"]) ./ repeat(reshape(mapslices(sum, (start_data["evfp"]), dims=1)[1, :, :], [1, size(start_data["evfp"])[2], size(start_data["evfp"])[3]]...), inner=[size(start_data["evfp"])[1], 1, 1])
    calibrate_start["σ_vif"] = NamedArray(mapslices(sum, start_data["vdfp"] .+ start_data["vmfp"], dims=[1])[1, :, :] ./ (mapslices(sum, start_data["vdfp"] .+ start_data["vmfp"], dims=[1])[1, :, :] .+ mapslices(sum, start_data["evfp"], dims=[1])[1, :, :]), (sets["acts"], sets["reg"]))


    # CAL-III
    fixed["α_qfdqfm"] .= false
    fixed["α_qfa"] .= false
    fixed["σ_vf"] = NamedArray(trues(size(data["σ_vf"])), names(data["σ_vf"]))
    fixed["σ_vf"][1, :, :] .= false
    fixed["σ_vdf"] = NamedArray(trues(size(data["σ_vdf"])), names(data["σ_vdf"]))
    fixed["ϵ_qfdqfm"] = NamedArray(trues(size(data["ϵ_qfdqfm"])), names(data["ϵ_qfdqfm"]))
    fixed["ϵ_qfa"] = NamedArray(trues(size(data["ϵ_qfa"])), names(data["ϵ_qfa"]))

    calibrate_start["σ_vf"] = (start_data["vdfp"] .+ start_data["vmfp"]) ./ repeat(reshape(mapslices(sum, (start_data["vdfp"] .+ start_data["vmfp"]), dims=1)[1, :, :], [1, size(start_data["vmfp"])[2], size(start_data["vmfp"])[3]]...), inner=[size(start_data["vmfp"])[1], 1, 1])
    calibrate_start["σ_vdf"] = (start_data["vdfp"]) ./ (start_data["vdfp"] .+ start_data["vmfp"])


    # CAL-IV
    fixed["α_qpdqpm"] .= false
    fixed["β_qpa"] .= false
    fixed["ϵ_qpdqpm"] = NamedArray(trues(size(data["ϵ_qpdqpm"])), names(data["ϵ_qpdqpm"]))
    fixed["σ_vp"] = NamedArray(trues(size(data["σ_vp"])), names(data["σ_vp"]))
    fixed["σ_vp"][1, :] .= false
    fixed["σ_vdp"] = NamedArray(trues(size(data["σ_vdp"])), names(data["σ_vdp"]))
    fixed["u"] = NamedArray(trues(size(data["u"])), names(data["u"])...)

    calibrate_start["σ_vp"] = (start_data["vdpp"] .+ start_data["vmpp"]) ./ repeat(reshape(mapslices(sum, (start_data["vdpp"] .+ start_data["vmpp"]), dims=1)[1, :], [1, size(start_data["vmpp"])[2]]...), inner=[size(start_data["vmpp"])[1], 1])
    calibrate_start["σ_vdp"] = (start_data["vdpp"]) ./ (start_data["vdpp"] .+ start_data["vmpp"])

    # CAL-V
    fixed["α_qgdqgm"] .= false
    fixed["α_qga"] .= false
    fixed["σ_vg"] = NamedArray(trues(size(data["σ_vg"])), names(data["σ_vg"]))
    fixed["σ_vg"][1, :] .= false
    fixed["σ_vdg"] = NamedArray(trues(size(data["σ_vdg"])), names(data["σ_vdg"]))
    fixed["ϵ_qgdqgm"] = NamedArray(trues(size(data["ϵ_qgdqgm"])), names(data["ϵ_qgdqgm"]))
    fixed["ϵ_qga"] = NamedArray(trues(size(data["ϵ_qga"])), names(data["ϵ_qga"])...)

    calibrate_start["σ_vg"] = (start_data["vdgp"] .+ start_data["vmgp"]) ./ repeat(reshape(mapslices(sum, (start_data["vdgp"] .+ start_data["vmgp"]), dims=1)[1, :], [1, size(start_data["vmgp"])[2]]...), inner=[size(start_data["vmgp"])[1], 1])
    calibrate_start["σ_vdg"] = (start_data["vdgp"]) ./ (start_data["vdgp"] .+ start_data["vmgp"])


    # CAL-VI
    fixed["α_qidqim"] .= false
    fixed["α_qia"] .= false
    fixed["σ_vi"] = NamedArray(trues(size(data["σ_vi"])), names(data["σ_vi"]))
    fixed["σ_vi"][1, :] .= false
    fixed["σ_vdi"] = NamedArray(trues(size(data["σ_vdi"])), names(data["σ_vdi"]))
    fixed["ϵ_qidqim"] = NamedArray(trues(size(data["ϵ_qidqim"])), names(data["ϵ_qidqim"]))
    fixed["ϵ_qia"] = NamedArray(trues(size(data["ϵ_qia"])), names(data["ϵ_qia"])...)

    calibrate_start["σ_vi"] = (start_data["vdip"] .+ start_data["vmip"]) ./ repeat(reshape(mapslices(sum, (start_data["vdip"] .+ start_data["vmip"]), dims=1)[1, :], [1, size(start_data["vmip"])[2]]...), inner=[size(start_data["vmgp"])[1], 1])
    calibrate_start["σ_vdi"] = (start_data["vdip"]) ./ (start_data["vdip"] .+ start_data["vmip"])

    # CAL-VII
    fixed["α_qtmfsd"] .= false
    fixed["σ_vtwr"] = NamedArray(trues(size(data["σ_vtwr"])), names(data["σ_vtwr"]))
    calibrate_start["σ_vtwr"] = (start_data["vtwr"]) ./ repeat(reshape(start_data["vcif"], [1, size(start_data["vcif"])...]...), inner=[size(start_data["vtwr"])[1], 1, 1, 1])

    # CAL-VIII
    fixed["α_qinv"] .= false
    fixed["σ_qinv"] = NamedArray(trues(size(data["σ_qinv"])), names(data["σ_qinv"])...)
    fixed["σ_qinv"][1] = false
    fixed = merge(fixed, Dict("ϵ_qinv" => true))
    calibrate_start["σ_qinv"] = (start_data["psave"] .* start_data["qsave"]) / sum(start_data["psave"] .* start_data["qsave"])


    # CAL-IX
    fixed["ρ"] = NamedArray(falses(size(data["kb"])), names(data["kb"])...)
    fixed["σ_ρ"] = NamedArray(trues(size(data["kb"])), names(data["kb"])...)
    calibrate_start["σ_ρ"] = mapslices(sum, start_data["evos"], dims=[1, 2])[1, 1, :] ./ start_data["vkb"]

    fixed["ppa"] .= false
    fixed["vdpp"] = NamedArray(falses(size(data["vdpp"])), names(data["vdpp"]))
    fixed["vdpp"][1, 1] = true

    calibrate_start["vdpp"] = start_data["vdpp"]

    #mc = model_container_struct(JuMP.Model(Ipopt.Optimizer), calibrate_start, parameters, sets, fixed)
    #build_model!(mc)
    #run_model!(model_container=mc, max_iter=max_iter, constr_viol_tol=constr_viol_tol)

    return (fixed_calibration=deepcopy(fixed), data_calibration=deepcopy(calibrate_start))
end