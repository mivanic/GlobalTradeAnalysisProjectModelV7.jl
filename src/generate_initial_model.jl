function generate_initial_model(; hSets, hData, hParameters)
    (; sets=sets, parameters, data, fixed) = generate_starting_values(; hSets, hData, hParameters)
    mc = model_container(JuMP.Model(Ipopt.Optimizer), data, parameters, sets, fixed)
    build_model!(mc)
    return mc
end