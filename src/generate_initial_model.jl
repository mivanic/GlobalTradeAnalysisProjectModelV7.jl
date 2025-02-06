"""
    generate_initial_model(; hSets, hData, hParameters)
"""
function generate_initial_model(; hSets, hData, hParameters)
    (; sets, parameters, data, fixed, lower, upper) = generate_starting_values(; hSets, hData, hParameters)
    mc = model_container_struct(JuMP.Model(Ipopt.Optimizer), data, parameters, sets, fixed, lower, upper)
    build_model!(mc)
    initialize_model!(mc)
    return mc
end