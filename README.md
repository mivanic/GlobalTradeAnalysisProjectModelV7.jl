# GTAPv7

[![Build Status](https://github.com/mivanic/GTAPv7.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/mivanic/GTAPv7.jl/actions/workflows/CI.yml?query=branch%3Amaster)

## Purpose of the package

The purpose of this package is to allow CGE modelers to run GTAP v7 in Julia, from start to finish, i.e. from data aggregation to results.

To run the GTAPv7 model, you need to follow the following steps:

## Prerequisites

- Get GTAP sets, data and parameters prepared for version 7 of the GTAP model
    - The sets, data and parameters need to be provided as a dictionary with the correct keys matching the headers of the model file for GEMPACK in lower case. For example, the data should contain key "vfob" with a three dimensional array with the value of bilateral trade measured FOB
    - If you are interested in actual GTAP data, you may obtain a free dataset from [the GTAP Center's website](https://www.gtap.agecon.purdue.edu/)
    - To turn a Fortran-style HAR file into a Julia dictionary, you may use package [HeaderArrayFile](https://github.com/mivanic/HeaderArrayFile.jl)
    - You may either aggregate your data outside module `GTAPv7`, e.g., by using GTAPAgg or FlexAgg programs also distributed by the GTAP Center, or you can use function `aggregate_data` in model `GTAPv7` as explained below
- A tiny aggregation (four regions, four commodities, four factors) of the GTAP database is provided with the package for testing using function get_sample_data()

## Aggregating the data

- To aggregate GTAP data in module `GTAPv7`, you can run function `aggregate_data(; hData, hParameters, hSets, comMap, regMap, endMap)` with the following arguments:
    - `hData`: a Dict object with the keys appropriate for the GTAPv7 model data (e.g., "vfob", "evos", etc.)
    - `hParameters`: a Dict object with the keys appropriate for the GTAPv7 model parameters (e.g., "esbq", etc.)
    - `comMap`: a Named Vector which maps original commodities to the new ones 
    - `regMap`: a Named Vector which maps original regions to the new ones
    - `endMap`: a Named Vector which maps original endowments to the new ones
        

```
# Use packages HeaderArrayFile and NamedArrays to process the initial data
using HeaderArrayFile, NamedArrays

# Load the disaggregated data (e.g., from FlexAgg)
parameters = HeaderArrayFile.readHar("./gsdfpar.har")
data = HeaderArrayFile.readHar("./gsdfdat.har")
sets = HeaderArrayFile.readHar("./gsdfset.har")

# Aggregate the data

## Prepare the mapping vectors
## You can start by reading the regions, commodities and endowments from the disaggregated data and modifying them
regMap = NamedArray(sets["reg"], sets["reg"])
regMap[1:27] .= "oceania+asia"
regMap[28:55] .= "americas"
regMap[56:98] .= "eu+oeurope"
regMap[99:160] .= "mena+africa"


comMap = NamedArray(copy(sets["comm"]), copy(sets["comm"]))
comMap[1:12] .= "ag"
comMap[13:18] .= "extract"
comMap[19:45] .= "manuf"
comMap[46:65] .= "svces"

endMap = NamedArray(copy(sets["endw"]), copy(sets["endw"]))
endMap[:] .= "other"
endMap[1:1] .= "land"
endMap[[2,5]] .= "labor"
endMap[[3,4,6]] .= "labor"
endMap["capital"] = "capital"

using GeneralEquilibrium.ModelLibrary.GTAPv7

# Do the aggregation
(; hData, hParameters, hSets) = GTAPv7.aggregate_data(hData=data, hParameters=parameters, hSets=sets, comMap=comMap, regMap=regMap, endMap=endMap)

```

# Alternatively, get a sample dataset

```

(; hData, hParameters, hSets) = GTAPv7.get_sample_data()

```

# Generating starting values

```
# Starting data and parameters
(; sets, parameters, data, fixed) = GTAPv7.generate_starting_values(hSets=hSets, hData=hData, hParameters=hParameters)
start_data = deepcopy(data)
```

# Solve the model with the starting (uncalibrated) data  values

```
(; data) = GTAPv7.model(sets=sets, data=start_data, parameters=parameters, fixed=fixed, max_iter=30, constr_viol_tol = 1e-8)
calibrated_data = deepcopy(data)
```


# Calibrate the data and parameters

```
(; data, parameters) = GTAPv7.calibrate(start_data = start_data, data=data, sets=sets,  start_parameters=parameters, fixed=fixed)

calibrated_data = copy(data)
```


# Running baseline scenario (the world before the shock)

```
(; data) = GTAPv7.model(sets=sets, data=calibrated_data, parameters=parameters,  fixed=fixed, calibrate=false, max_iter = 20)

# Let's save the state of the world before the simulation
data0 = deepcopy(data)
```

# Running the scenario (the world after the shock)

```
# Set the tariff on crops from ssafrica to eu to 1.2 (20 percent)
data["tms"]["crops", "mena+africa", "eu"] = 1.2

# Run the model
(; data, calibrated_parameters) = GTAPv7.model(sets=sets, data=data, parameters=parameters, calibrated_parameters=calibrated_calibrated_parameters, fixed=fixed, hData=hData, calibrate=false, max_iter=20)

# Save the world after the simulation
data1 = deepcopy(data)
```

# Analyzing the results

```
# Show the change in exports (percent)
((data1["qxs"]./data0["qxs"])[:, :, "eu+oeurope"] .- 1) .* 100

```




