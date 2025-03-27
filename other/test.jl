include("./../src/GlobalTradeAnalysisProjectModelV7.jl")
using Main.GlobalTradeAnalysisProjectModelV7
using HeaderArrayFile, NamedArrays
#using BenchmarkTools


parameters = HeaderArrayFile.readHar("./../GTAPDatabase/gsdfpar.har")
data = HeaderArrayFile.readHar("./../GTAPDatabase/gsdfdat.har")
sets = HeaderArrayFile.readHar("./../GTAPDatabase/gsdfset.har")

#data["vtwr"][:,"ctl","mli","tgo"]
#data["vcif"]["ctl","mli","tgo"]


regMap = NamedArray(deepcopy(sets["reg"]), sets["reg"])
regMap[1:3] .= "Oceania"
regMap[4] = "China"
regMap[5:27] .= "Rest_Asia"
regMap[28] = "Canada"
regMap[29] = "USA"
regMap[30] = "Mexico"
regMap[31] = "RestofWorld"; regMap[94:101].="RestofWorld"; regMap[160] = "RestofWorld"
regMap[32:55] .= "LatinAmer"
regMap[56:82] .= "EU"
regMap[83:89] .= "Rest_Eur"; regMap[91:93] .= "Rest_Eur"
regMap[90] = "Russia"
regMap[102:121] .= "MENA"
regMap[122:159] .= "SSA"

# regMap[4] = "Rest_Asia"
# regMap[56:82] .= "Rest_Eur"
# regMap[90] = "Rest_Eur" ## SUSPICOUS!
# regMap[102:121] .= "RestofWorld"
# regMap[122:159] .= "RestofWorld"
# regMap[28:29] .= "RestofWorld"
# regMap[1:4] .= "RestofWorld"
# regMap[regMap.=="Rest_Asia"] .= "RestofWorld"
# regMap[regMap.=="Rest_Eur"] .= "RestofWorld"
# regMap[regMap.=="Mexico"] .= "RestofWorld"

#regMap["rus"] = "Russia"

# regMap[1:3] .= "oceania"
# regMap[4:10] .= "easia"
# regMap[11:20] .= "seasia"
# regMap[21:27] .= "sasia"
# regMap[28:31] .= "namerica"
# regMap[32:42] .= "samerica"
# regMap[43:55] .= "carib"
# regMap[56:82] .= "eu"
# regMap[83:98] .= "oeurope"
# regMap[99:116] .= "wasia"
# regMap[117:121] .= "mena"
# regMap[122:133] .= "wafrica"
# regMap[134:140] .= "scafrica"
# regMap[141:154] .= "eafrica"
# regMap[155:160] .= "safrica"

# regMap[1:3] .= "oceania"
# regMap[4:10] .= "easia"
# regMap[11:20] .= "seasia"
# regMap[21:27] .= "sasia"
# regMap[28:31] .= "namerica"
# regMap[32:42] .= "samerica"
# regMap[43:55] .= "carib"
# regMap[56:82] .= "eu"
# regMap[83:98] .= "oeurope"
# regMap[99:116] .= "wasia"
# regMap[117:121] .= "mena"
# regMap[122:160] .= "wafrica"

comMap = NamedArray(deepcopy(sets["comm"]), sets["comm"])
comMap[1:8] .= "GrainsCrops"
comMap[9:12] .= "MeatLstk"
comMap[13:14] .= "Extraction"
comMap[15] = "Coal"
comMap[16] = "Oil"
comMap[17] = "Gas"
comMap[18] = "Extraction"
comMap[19:20] .= "MeatLstk"
comMap[21:22] .= "ProcFood"
comMap[23] = "GrainsCrops"
comMap[24:26] .= "ProcFood"
comMap[27:28] .= "TextWapp"
comMap[29:31] .= "LightMnfc"
comMap[32:38] .= "HeavyMnfc"
comMap[39] = "LightMnfc"
comMap[40:42] .= "HeavyMnfc"
comMap[43] = "Vehicles"
comMap[44:45] .= "LightMnfc"
comMap[46:49] .= "Util_Cons"
comMap[50:56] .= "TransComm"
comMap[57:65] .= "OthServices"

# comMap[19:65] .= "Rest"
# comMap[1:12] .= "Ag"
# comMap[[15,16]] .= "Extraction"

# comMap[1:8] .= "crops"
# comMap[9:12] .= "animals"
# comMap[13:18] .= "extract"
# #comMap[15]="Extraction"
# #comMap[16]="Extraction"
# #comMap[17]="Extraction"
# comMap[19:26] .= "food"
# comMap[27:45] .= "manuf"
# comMap[46:65] .= "svces"

endMap = NamedArray(deepcopy(sets["endw"]), sets["endw"])
endMap[1:1] .= "Land"
endMap[[2, 5]] .= "SkLab"
endMap[[3, 4, 6]] .= "UnSkLab"
endMap["capital"] = "capital"
endMap[8:8] .= "NatlRes"

# Do the aggregation
(; hData, hParameters, hSets) = aggregate_data_legacy(; hData=data, hParameters=parameters, hSets=sets, comMap=comMap, regMap=regMap, endMap=endMap)

# Clean the data
#(; hData, hParameters, hSets) = clean_data(; hData=hData, hParameters=hParameters, hSets=hSets)

orig_hData = deepcopy(hData)
#hParameters["esbm"][hParameters["esbm"].>10].=10 #34.4
#hParameters["esbd"][hParameters["esbd"].>5].=5 #17.2

#hParameters["esbm"][hParameters["esbm"].>5] .= 5
#hParameters["esbd"]["Gas",:] .= hParameters["esbd"]["Gas",:] ./ 5

# # Clean the data: trade below 1e-6 makes little sense
# valid_trade = hData["vxsb"].>1e-1
# for k ∈ ["vcif","vmsb","vfob","vxsb"]
#     hData[k][.!valid_trade].=0
# end


# for k ∈ ["vtwr"]
#     for m ∈ hSets["marg"]
#         for c ∈ hSets["comm"]
#             for s ∈ hSets["reg"]
#                 for d ∈ hSets["reg"]
#                     if !valid_trade[c, s, d]
#                         hData[k][m,c,s,d] .= 0
#                     end
#                 end
#             end
#         end
#     end
# end
# valid_margins = hData["vtwr"].>1e-8
# for k ∈ ["vtwr"]
#     hData[k][.!valid_margins].=0
# end

# valid_inputs = hData["vdfp"] .+ hData["vmfp"] .>10
# # # for c ∈ hSets["comm"]
# # #     for a ∈ hSets["acts"]
# # #         for r ∈ hSets["reg"]
# # #             if (hData["vdfp"][c,a,r] + hData["vmfp"][c,a,r])/sum(hData["vdfp"][:,a,r] .+ hData["vmfp"][:,a,r]) < 1e-2
# # #                 valid_inputs[c,a,r]  = false
# # #             end
# # #         end
# # #     end
# # # end

# hData["vdfp"][.!valid_inputs] .=0
# hData["vmfp"][.!valid_inputs] .=0
# hData["vdfb"][.!valid_inputs] .=0
# hData["vmfb"][.!valid_inputs] .=0

# #minimum(hData["vdfp"] .+ hData["vmfp"])
# sum(.!valid_inputs)


# valid_gconsumption = hData["vdgp"] .+ hData["vmgp"] .>10
# # # for c ∈ hSets["comm"]
# # #     for a ∈ hSets["acts"]
# # #         for r ∈ hSets["reg"]
# # #             if (hData["vdfp"][c,a,r] + hData["vmfp"][c,a,r])/sum(hData["vdfp"][:,a,r] .+ hData["vmfp"][:,a,r]) < 1e-2
# # #                 valid_inputs[c,a,r]  = false
# # #             end
# # #         end
# # #     end
# # # end

# hData["vdgp"][.!valid_gconsumption] .=0
# hData["vmgp"][.!valid_gconsumption] .=0
# hData["vdgb"][.!valid_gconsumption] .=0
# hData["vmgb"][.!valid_gconsumption] .=0

# #minimum(hData["vdfp"] .+ hData["vmfp"])
# sum(.!valid_gconsumption)


#hParameters["esbm"][hParameters["esbm"].>10].=10
#hParameters["esbd"][hParameters["esbm"].>5].=5
mc = generate_initial_model(hSets=hSets, hData=hData, hParameters=hParameters)


start_data = deepcopy(mc.data)


# using DataFrames, ClipData, JuMP
# cliptable(DataFrame(name = name.(all_variables(mc.model)), start_value = start_value.(all_variables(mc.model))))

# run_model!(mc; max_iter = 40)

# cliptable(DataFrame(name = name.(all_variables(mc.model)), value = value.(all_variables(mc.model))))

# hData["vtwr"]["svces","extract",:,:]
# hData["vcif"]["extract",:,"wafrica"]
# hData["vfob"]["extract",:,"wafrica"]
# mc.data["qxs"]["extract",:,"wafrica"]

# mc.data["α_qxs"]["extract",:,"wafrica"]
# mc.data["pcif"]["extract",:,"wafrica"]
# hParameters["esbm"]

# value.(mc.model[:qxs]["Coal",:,"LatinAmer"])
# value.(mc.model[:pcif]["Coal",:,"LatinAmer"])
# mc.data["qtmfsd"][:,"Coal",:,"LatinAmer"]
# mc.data["pt"]
# value.(mc.model[:pcif]["Coal",:,"LatinAmer"])
# mc.data["σ_qxs"]["Coal",:,"LatinAmer"]
# mc.data["qxs"]["Coal",:,"LatinAmer"] .* mc.data["pcif"]["Coal",:,"LatinAmer"]
# start_data["qxs"]["Coal",:,"LatinAmer"] .* start_data["pcif"]["Coal",:,"LatinAmer"]

# start_data["qms"].*start_data["pms"]
# mc.data["qms"].*mc.data["pms"]

# (start_data["qfm"].*start_data["pfm"])[:,:,"LatinAmer"]
# (start_data["qfd"].*start_data["pfd"])[:,:,"LatinAmer"]
# (start_data["qpm"].*start_data["ppm"])
# (start_data["qgm"].*start_data["pgm"])


# hData["vmfp"][:,:,"LatinAmer"]
# hData["vmfb"][:,:,"LatinAmer"]

# hData["vmfp"]["Coal",:,"LatinAmer"]
# hData["vmpp"]["Coal","LatinAmer"]
# hData["vmgp"]["Coal","LatinAmer"]
# hData["vmip"]["Coal","LatinAmer"]
# hData["vmfb"]["Coal",:,"LatinAmer"]
# hData["vmpb"]["Coal","LatinAmer"]
# hData["vmgb"]["Coal","LatinAmer"]
# hData["vmib"]["Coal","LatinAmer"]
# hData["vcif"]["Coal",:,"LatinAmer"]
# hData["vfob"]["Coal",:,"LatinAmer"]

# value.(mc.model[:up] .^ mc.model[:σyp] .* mc.model[:ug] .^ mc.model[:σyg] .* mc.model[:us] .^ (1 .- mc.model[:σyp] .- mc.model[:σyg]))

(;fixed_calibration, data_calibration)=generate_calibration_inputs(mc, start_data)

fixed_default = deepcopy(mc.fixed)

mc.data = data_calibration
mc.fixed = fixed_calibration

run_model!(mc)
#run_model!(mc)
#run_model!(mc)

#cliptable(DataFrame(name = name.(all_variables(mc.model)), value = value.(all_variables(mc.model))))

calibrated_data = deepcopy(mc.data)

mc.fixed = deepcopy(fixed_default)


rebuild_model!(mc)

@time run_model!(mc)


# using JuMP, Ipopt
# mc.model = JuMP.Model(Ipopt.Optimizer)
# build_model!(mc; calibration = false)


# # Delete unneded variables/equations
# using JuMP 
# delete.(mc.model, mc.model[:sf_α_qxs])
# delete.(mc.model, mc.model[:ϵ_qxs])
# delete.(mc.model,mc.model[:sf_α_qfe])
# delete.(mc.model,mc.model[:ϵ_qfe])
# delete.(mc.model,mc.model[:sf_α_qes2])
# delete.(mc.model,mc.model[:ϵ_qes2])
# delete.(mc.model,mc.model[:sf_α_qfdqfm])
# delete.(mc.model,mc.model[:ϵ_qfdqfm])

# #length(mc.model[:e_σ_qxs])
# #length(mc.model[:σ_qxs])
# delete.(mc.model, mc.model[:e_σ_qxs])
# delete.(mc.model, mc.model[:σ_qxs])
# delete.(mc.model, mc.model[:e_σ_vtwr])
# delete.(mc.model, Array(mc.model[:σ_vtwr])[Array(is_valid.(mc.model, mc.model[:σ_vtwr]))])

# delete.(mc.model, mc.model[:e_σ_vff])
# delete.(mc.model, Array(mc.model[:σ_vff])[Array(is_valid.(mc.model, mc.model[:σ_vff]))])

# delete.(mc.model, mc.model[:e_σ_vf])
# delete.(mc.model, Array(mc.model[:σ_vf])[Array(is_valid.(mc.model, mc.model[:σ_vf]))])

# delete.(mc.model, mc.model[:e_σ_vdf])
# delete.(mc.model, Array(mc.model[:σ_vdf])[Array(is_valid.(mc.model, mc.model[:σ_vdf]))])


# # 27627
# # 25777 - 1350 = 24427

# @time run_model!(mc)



#mc.data["tms"]["GrainsCrops","RestofWorld","LatinAmer"]=2
mc.data["tms"]["GrainsCrops","MENA","EU"]=2

@time run_model!(mc)

round.((mc.data["qxs"]["crops",:,:]./calibrated_data["qxs"]["crops",:,:] .- 1) .*100, digits = 1)

mc.data["walras_dem"]
mc.data["walras_sup"]

unique(regMap)
unique(comMap)
