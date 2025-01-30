function initialize_model!(mc)
    # First set the starting values and report if anything is missing
    for k ∈ names(object_dictionary(mc.model))
        if mc.model[k] isa VariableRef
            if String(k) ∈ names(mc.data)
                set_start_value(mc.model[k], mc.data[String(k)])
            else
                printstyled("Values for $k are not provided in the data\n", color=:yellow)
            end
        elseif mc.model[k] isa ConstraintRef
        elseif typeof(Array(mc.model[k])) ∈ [Array{VariableRef,1}, Array{VariableRef,2}, Array{VariableRef,3}, Array{VariableRef,4}, Array{VariableRef,5}]
            if String(k) ∈ names(mc.data)
                set_start_value.(Array(mc.model[k])[Array(is_valid.(mc.model, mc.model[k]))], mc.data[String(k)][Array(is_valid.(mc.model, mc.model[k]))])
            else
                printstyled("Values for $k are not provided in the data\n", color=:yellow)
            end
        else
        end
    end
    # Second unfix all variables 
    unfix.(JuMP.all_variables(mc.model)[is_fixed.(JuMP.all_variables(mc.model))])
    # Third fix variables that need to be fixed
    for fv ∈ keys(mc.fixed)
        if size(mc.fixed[fv]) == ()
            if mc.fixed[fv] && is_valid(mc.model, mc.model[Symbol(fv)])
                if !isnan(mc.data[fv])
                    fix(mc.model[Symbol(fv)], mc.data[fv]; force=true)
                else
                    fix(mc.model[Symbol(fv)], 0; force=true)
                end
            end
        else
            for fvi ∈ CartesianIndices(mc.fixed[fv])
                if mc.fixed[fv][fvi] && is_valid(mc.model, mc.model[Symbol(fv)][fvi])
                    if !isnan( mc.data[fv][fvi]) 
                        fix(mc.model[Symbol(fv)][fvi], mc.data[fv][fvi]; force=true)
                    else
                        fix(mc.model[Symbol(fv)][fvi], 0; force=true)
                    end
                end
            end
        end
    end
    return nothing
end