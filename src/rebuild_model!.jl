"""
    rebuild_model(mc; calibration = false)

Rebuilds the model, typically to create a slimmer model without calibration equations

### Args:
- `mc`: model container
- `calibration`: include calibration equations

## Returns:

An updated model container

### Example:

``julia
rebuild_model!(mc)
``

"""
function rebuild_model!(mc; model=JuMP.Model(Ipopt.Optimizer), calibration = false)
    mc.model = model
    build_model!(mc; calibration = calibration)
    return mc
end