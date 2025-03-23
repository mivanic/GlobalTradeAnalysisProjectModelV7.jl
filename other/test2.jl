include("./src/GlobalTradeAnalysisProjectModelV7.jl")
using Main.GlobalTradeAnalysisProjectModelV7
using HeaderArrayFile, NamedArrays
using BenchmarkTools


parameters = HeaderArrayFile.readHar("./../GtapInJump/Replication/Data/gsdfpar.har")
data = HeaderArrayFile.readHar("./../GtapInJump/Replication/Data/gsdfdat.har")
sets = HeaderArrayFile.readHar("./../GtapInJump/Replication/Data/gsdfset.har")

regMap = NamedArray(deepcopy(sets["reg"]), sets["reg"])
regMap .= "Rest of the World"
regMap["usa"]="usa"
regMap["can"]="can"
regMap["chn"]="chn"
regMap["mex"]="mex"

comMap = NamedArray(deepcopy(sets["comm"]), sets["comm"])
comMap[1:8] .= "Crops"
comMap[27:45] .= "Manufacturing"
comMap[46:65] .= "Services"
comMap[["coa","oil","gas","oxt"]] .= "Energy & Exctraction"
#comMap[["coa","oil","gas"]] .= "Energy"


comMap[["pdr","pcr"]].="Rice"
comMap[["ctl","cmt"]].="Cattle & meat"
comMap[["oap","omt","wol"]].="Other animals, meat & wool"
comMap[["c_b","sgr"]].="Sugar"
comMap[["rmk","mil"]].="Milk"
comMap[["tex","wap","lea"]].="Textiles & clothes"

unique(comMap)

endMap = NamedArray(deepcopy(sets["endw"]), sets["endw"])
endMap[1:1] .= "Land"
endMap[[2, 5]] .= "SkLab"
endMap[[3, 4, 6]] .= "UnSkLab"
endMap["capital"] = "capital"
endMap[8:8] .= "NatlRes"

# Do the aggregation
(; hData, hParameters, hSets) = aggregate_data_legacy(hData=data, hParameters=parameters, hSets=sets, comMap=comMap, regMap=regMap, endMap=endMap)


# Clean the data: trade below 1e-6 makes little sense
valid_trade = hData["vxsb"].>1e-6
for k ∈ ["vcif","vmsb","vfob","vxsb"]
    hData[k][.!valid_trade].=0
end

for k ∈ ["vtwr"]
    for m ∈ hSets["marg"]
        for c ∈ hSets["comm"]
            for s ∈ hSets["reg"]
                for d ∈ hSets["reg"]
                    if !valid_trade[c, s, d]
                        hData[k][m,c,s,d] .= 0
                    end
                end
            end
        end
    end
end
valid_margins = hData["vtwr"].>1e-6
for k ∈ ["vtwr"]
    hData[k][.!valid_margins].=0
end

valid_inputs = hData["vdfp"] .+ hData["vmfp"] .>1e-6
for c ∈ hSets["comm"]
    for a ∈ hSets["acts"]
        for r ∈ hSets["reg"]
            if (hData["vdfp"][c,a,r] + hData["vmfp"][c,a,r])/sum(hData["vdfp"][:,a,r] .+ hData["vmfp"][:,a,r]) < 1e-6
                valid_inputs[c,a,r]  = false
            end
        end
    end
end

hData["vdfp"][.!valid_inputs] .=0
hData["vmfp"][.!valid_inputs] .=0


mc = generate_initial_model(hSets=hSets, hData=hData, hParameters=hParameters)

start_data = deepcopy(mc.data)
run_model!(mc; max_iter = 20)


using JuMP
mc.data["σ_vtwr"]

using JuMP, ClipData, DataFrames
cliptable(DataFrame(name = name.(all_constraints(mc.model; include_variable_in_set_constraints = false )),value = value.(all_constraints(mc.model; include_variable_in_set_constraints = false ))))

value.(mc.model[:e_qgdqgm]["fsh","mex"])
mc.model[:e_qgdqgm]["fsh","mex"][2]

value.(mc.model[:qgm]["fsh","mex"])
value.(mc.model[:α_qgdqgm][:,"fsh","mex"])
value.(mc.model[:γ_qgdqgm]["fsh","mex"])

value((mc.model[:qgm]["fsh","mex"]) )
value(
    #log(mc.model[:qgm]["fsh","mex"]) 
    - log((mc.model[:qga]["fsh","mex"] / mc.model[:γ_qgdqgm]["fsh","mex"]) * ((((mc.model[:α_qgdqgm]["imp","fsh","mex"]*mc.model[:γ_qgdqgm]["fsh","mex"]) * ((1.0 / mc.model[:γ_qgdqgm]["fsh","mex"]) * (((+(0.0) + ((mc.model[:α_qgdqgm]["dom","fsh","mex"] ^ 1.2499998807907104) * (mc.model[:pgd]["fsh","mex"] ^ -0.24999988079071045))) + ((mc.model[:α_qgdqgm]["imp","fsh","mex"] ^ 1.2499998807907104) * (mc.model[:pgm]["fsh","mex"] ^ -0.24999988079071045))) ^ -4.000001907349542))) / mc.model[:pgm]["fsh","mex"]) ^ 1.2499998807907104)) )

a="Services"
r="Rest of the World"

log(σ_vif[a, r]) + log(pva[a, r] * qva[a, r] + pint[a, r] * qint[a, r]) == log(pint[a, r] * qint[a, r])
value(log(mc.model[:σ_vif][a, r]) + log(mc.model[:pva][a, r] * mc.model[:qva][a, r] + mc.model[:pint][a, r] * mc.model[:qint][a, r]))
value(log(mc.model[:pint][a, r] * mc.model[:qint][a, r]))

value(mc.model[:e_σ_vif][a,r])

cliptable(DataFrame(name = name.(all_variables(mc.model)),value = start_value.(all_variables(mc.model))))


run_model!(mc; max_iter = 0)
run_model!(mc; max_iter = 30)
run_model!(mc; max_iter = 30)
run_model!(mc; max_iter = 30)
run_model!(mc; max_iter = 30)
run_model!(mc; max_iter = 30)
run_model!(mc; max_iter = 30)

(;fixed_calibration, data_calibration)=generate_calibration_inputs(mc, start_data)

fixed_default = deepcopy(mc.fixed)

mc.data = data_calibration
mc.fixed = fixed_calibration

run_model!(mc)

calibrated_data = deepcopy(mc.data)

mc.fixed = deepcopy(fixed_default)

run_model!(mc)



using JuMP, ClipData
cliparray(name.(all_constraints(mc.model; include_variable_in_set_constraints = false )))

unique(comMap)
unique(regMap)
