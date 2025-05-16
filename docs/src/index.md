# Introduction

This package is intended to allow the user to run the standard GTAP model v 7 in levels and completely in Julia (JuMP, Ipopt).


Install:

```julia
using Pkg
Pkg.add("GlobalTradeAnalysisProjectModelV7")
```

Useful packages:

```julia
using Pkg
Pkg.add("HeaderArrayFile")
Pkg.add("NamedArrays")
```

## Simple example

This example will use a sample data set to do a simple tariff simulation.

```julia
# Load packages
using GlobalTradeAnalysisProjectModelV7, HeaderArrayFile, NamedArrays

# Get the sample data
(hData, hParameters, hSets) = get_sample_data()

# Clean the data: trade below 1e-6 makes little sense
valid_trade = hData["vxsb"].>1e-4
for k ∈ ["vcif","vmsb","vfob","vxsb"]
    hData[k][.!valid_trade].=0
end

# Produce initial uncalibrated model using the GTAP data
mc = generate_initial_model(hSets=hSets, hData=hData, hParameters=hParameters)

# Keep the start data for calibration---the value flows are the correct ones
start_data = deepcopy(mc.data)

# Get the required inputs for calibration by providing the target values in start_data
(;fixed_calibration, data_calibration)=generate_calibration_inputs(mc, start_data)

# Keep the default closure (fixed) for later
fixed_default = deepcopy(mc.fixed)

# Load the calibration data and closure 
mc.data = deepcopy(data_calibration)
mc.fixed = deepcopy(fixed_calibration)

run_model!(mc)

# Save the calibrated data---this is the starting point for all simulation
calibrated_data = deepcopy(mc.data)

# Let's change the closure to the default (simulation) closure
mc.fixed = deepcopy(fixed_default)

# Start with the calibrated data
mc.data = deepcopy(calibrated_data)


### TARIFF SCENARIO
# Double the power of tariff
mc.data["tms"][["crops", "processed food"], ["mena", "sub-saharan africa"], "eu"] .= mc.data["tms"][["crops", "processed food"], ["mena", "sub-saharan africa"], "eu"] * 2

# Run the model
run_model!(mc)

## View some of the solutions:
# See change in exports to the eu
round.((mc.data["qxs"][:,:,"eu"] ./ calibrated_data["qxs"][:,:,"eu"] .-1) .* 100,digits = 2)

# See the percentage change in qpa, for example:
round.((mc.data["qpa"] ./ calibrated_data["qpa"] .-1) .* 100, digits = 2)

# Calculate EV
ev = calculate_expenditure(sets = mc.sets, data0=calibrated_data, data1=mc.data, parameters=mc.parameters)  .- calibrated_data["y"]
```

## Changes to the model


Because the model container contains an actual JuMP model, it is possible to make changes to the model without rewriting it. 


### Example

Let's modify the assumption of the standard model that land supply is exogenous by replacing this assumption with a constant elasticity supply equation:

``\mbox{qe}_{\mbox{land}}=\alpha_{\mbox{land}}\times\left(\frac{\mbox{pe}_{\mbox{land}}}{\mbox{ppriv}}\right)^{\sigma_{\mbox{land}}}``

Add new variables:

```julia
reg= mc.sets["reg"]
@variables(mc.model,
    begin
        land_scale[reg]
        land_elasticity[reg]
    end
)
```

Add the equation:

```julia 
@constraints(mc.model,
    begin
        e_qo_land[r=reg], log(mc.model[:qe]["land", r]) == 
            log(land_scale[r] * 
                    (mc.model[:pe]["land",r] / 
                       mc.model[:ppriv][r]) ^ land_elasticity[r])
    end
)
```

Calibrate the new  model (let me assume that the land supply elasticity is -1):

```julia
# Calibrate the new parameters
mc.data["land_elasticity"]= -1 .* NamedArray(ones(length(reg)), reg)
```

Calculate the calibrated value of the scaling factor:

```julia 
mc.data["land_scale"] = 
  mc.data["qe"]["land",:] ./ ((mc.data["pe"]["land",:] ./ 
                mc.data["ppriv"]).^mc.data["land_elasticity"])
```

Finally,  make `qe` of land endogenous, and the new parameters exogenous:

```julia
mc.fixed["land_elasticity"]=
  NamedArray(trues(size(mc.data["land_elasticity"])), reg)
mc.fixed["land_scale"]= NamedArray(trues(size(mc.data["land_scale"])), reg)
mc.fixed["qe"]["land",:].=false
```

You can run the modified model:

```julia
run_model!(mc)
```

