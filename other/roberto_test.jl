using GlobalTradeAnalysisProjectModelV7
using HeaderArrayFile, NamedArrays

# Get the raw GTAP v 11 data (NOT PROVIDED!)
parameters = HeaderArrayFile.readHar("./gsdfpar.har")
data = HeaderArrayFile.readHar("./gsdfdat.har")
sets = HeaderArrayFile.readHar("./gsdfset.har")

# Specify mappings

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


endMap = NamedArray(deepcopy(sets["endw"]), sets["endw"])
endMap[1:1] .= "Land"
endMap[[2, 5]] .= "SkLab"
endMap[[3, 4, 6]] .= "UnSkLab"
endMap["capital"] = "capital"
endMap[8:8] .= "NatlRes"


# Do the aggregation
(; hData, hParameters, hSets) = aggregate_data_legacy(; hData=data, hParameters=parameters, hSets=sets, comMap=comMap, regMap=regMap, endMap=endMap)

# Generate the model
mc = generate_initial_model(hSets=hSets, hData=hData, hParameters=hParameters)

# Save the starting data for calibration
start_data = deepcopy(mc.data)

# Create inputs for calibration
(;fixed_calibration, data_calibration)=generate_calibration_inputs(mc, start_data)

# Keep standard closure on the side for subsequent simulations
fixed_default = deepcopy(mc.fixed)

# Change the data with the calibration data and closure with calibration closure
mc.data = data_calibration
mc.fixed = fixed_calibration

# Run the model (calibration)
run_model!(mc)

# Save the calibrated data for subsequent simulations--this is the baseline
calibrated_data = deepcopy(mc.data)

# Change the closure back to the simulation closure
mc.fixed = deepcopy(fixed_default)

# Drop the equations that are not needed for solution
rebuild_model!(mc)

# Change the power of tariff to 2
mc.data["tms"]["GrainsCrops","MENA","EU"]=2

# Run the model
run_model!(mc)

# See some results
round.((mc.data["qxs"]["GrainsCrops",:,:]./calibrated_data["qxs"]["GrainsCrops",:,:] .- 1) .*100, digits = 3)
round.((mc.data["ppa"]./calibrated_data["ppa"] .- 1) .*100, digits = 3)

