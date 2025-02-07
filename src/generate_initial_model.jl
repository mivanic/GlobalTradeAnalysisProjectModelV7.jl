"""
    generate_initial_model(; hSets, hData, hParameters)

Produces an initial model container based on data, parameters and sets.

### Args:
- `hSet`: a dictionary of sets
- `hData`: a dictionary of data
- `hParameters`: a dictionary of parameters

## Returns:

A model container with the following elements:
- `model`: an empty Ipopt model
- `data`: a dictionary of data (variables in the model)
- `parameters`: a dictionary of parameter values
- `sets`: a dictionary of sets
- `fixed`: a dictionary of Boolean arrays defining the closure of the model (exogenous/fixed variables)
- `lower`: a dictionary of lower bounds for model variables
- `upper`: a dictionary of upper bounds for model variables

### Example:

``julia
generate_initial_model(; hSets, hData, hParameters)
``

"""
function generate_initial_model(; hSets, hData, hParameters)
    (; sets, parameters, data, fixed, lower, upper) = generate_starting_values(; hSets, hData, hParameters)
    mc = model_container_struct(JuMP.Model(Ipopt.Optimizer), data, parameters, sets, fixed, lower, upper)
    build_model!(mc)
    initialize_model!(mc)
    return mc
end