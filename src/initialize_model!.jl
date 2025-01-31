function initialize_model!(; model_container)

    free_variables = filter(f->!is_fixed(f), all_variables(model_container.model))
    set_start_value.(free_variables,1.01)

    # First set the starting values and report if anything is missing
    for k ∈ names(object_dictionary(model_container.model))
        if model_container.model[k] isa VariableRef
            if String(k) ∈ names(model_container.data)
                set_start_value(model_container.model[k], model_container.data[String(k)])
            else
                printstyled("Values for $k are not provided in the data\n", color=:yellow)
            end
        elseif model_container.model[k] isa ConstraintRef
        elseif typeof(Array(model_container.model[k])) ∈ [Array{VariableRef,1}, Array{VariableRef,2}, Array{VariableRef,3}, Array{VariableRef,4}, Array{VariableRef,5}]
            if String(k) ∈ names(model_container.data)
                set_start_value.(Array(model_container.model[k])[Array(is_valid.(model_container.model, model_container.model[k]))], model_container.data[String(k)][Array(is_valid.(model_container.model, model_container.model[k]))])
            else
                printstyled("Values for $k are not provided in the data\n", color=:yellow)
            end
        else
        end
    end
    # Second unfix all variables 
    unfix.(JuMP.all_variables(model_container.model)[is_fixed.(JuMP.all_variables(model_container.model))])
    # Third fix variables that need to be fixed
    for fv ∈ keys(model_container.fixed)
        if size(model_container.fixed[fv]) == ()
            if model_container.fixed[fv] && is_valid(model_container.model, model_container.model[Symbol(fv)])
                if !isnan(model_container.data[fv])
                    fix(model_container.model[Symbol(fv)], model_container.data[fv]; force=true)
                else
                    fix(model_container.model[Symbol(fv)], 0; force=true)
                end
            end
        else
            for fvi ∈ CartesianIndices(model_container.fixed[fv])
                if model_container.fixed[fv][fvi] && is_valid(model_container.model, model_container.model[Symbol(fv)][fvi])
                    if !isnan( model_container.data[fv][fvi]) 
                        fix(model_container.model[Symbol(fv)][fvi], model_container.data[fv][fvi]; force=true)
                    else
                        fix(model_container.model[Symbol(fv)][fvi], 0; force=true)
                    end
                end
            end
        end
    end
    return nothing
end