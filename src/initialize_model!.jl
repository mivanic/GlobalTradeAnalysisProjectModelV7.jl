function initialize_model!(mc)

    # Set starting values
    for k in keys(mc.data)
        if Symbol(k) ∈ keys(object_dictionary(mc.model))
            if mc.data[k] isa NamedArray
                set_start_value.(mc.model[Symbol(k)], Array(mc.data[k]))
                unfix.(Array(mc.model[Symbol(k)])[Array(is_fixed.(mc.model[Symbol(k)]))])
            else
                set_start_value.(mc.model[Symbol(k)], mc.data[k])
                if is_fixed.(mc.model[Symbol(k)])
                    unfix.(mc.model[Symbol(k)])
                end
            end
        end
    end
    
    # Fix fixed values and delete missing ones
    for fv ∈ keys(mc.fixed)
        if size(mc.fixed[fv]) == ()
            if mc.fixed[fv] && is_valid(mc.model, mc.model[Symbol(fv)])
                if isnan(mc.data[fv])
                    if is_valid(mc.model, mc.model[Symbol(fv)])
                        delete(mc.model, mc.model[Symbol(fv)])
                    end
                else
                    fix(mc.model[Symbol(fv)], mc.data[fv]; force=true)
                end
            end
        else
            for fvi ∈ CartesianIndices(mc.fixed[fv])
                if mc.fixed[fv][fvi] && is_valid(mc.model, mc.model[Symbol(fv)][fvi])
                    if isnan(mc.data[fv][fvi]) 
                        if is_valid(mc.model, mc.model[Symbol(fv)][fvi]) 
                            delete(mc.model, mc.model[Symbol(fv)][fvi])
                        end
                    else
                        fix(mc.model[Symbol(fv)][fvi], mc.data[fv][fvi]; force=true)
                    end
                end
            end
        end
    end
    return nothing
end