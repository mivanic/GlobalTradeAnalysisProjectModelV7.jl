function initialize_model!(mc)

    # Set starting values
    for k in keys(mc.data)
        if Symbol(k) ∈ keys(object_dictionary(mc.model))
            if mc.data[k] isa NamedArray
                set_start_value.(mc.model[Symbol(k)], Array(mc.data[k]))
                unfix.(mc.model[Symbol(k)], Array(mc.data[k]))
            else
                set_start_value.(mc.model[Symbol(k)], mc.data[k])
                unfix.(mc.model[Symbol(k)], mc.data[k])
            end
        end
    end

    # Fix fixed values and delete missing ones
    for fv ∈ keys(fixed)
        if size(fixed[fv]) == ()
            if fixed[fv] && is_valid(mc.model, mc.model[Symbol(fv)])
                if isnan(mc.data[fv])
                    delete(mc.model, mc.model[Symbol(fv)])
                else
                    fix(mc.model[Symbol(fv)], mc.data[fv]; force=true)
                end
            end
        else
            for fvi ∈ CartesianIndices(fixed[fv])
                if fixed[fv][fvi] && is_valid(mc.model, mc.model[Symbol(fv)][fvi])
                    if isnan(mc.data[fv][fvi])
                        delete(mc.model, mc.model[Symbol(fv)][fvi])
                    else
                        fix(mc.model[Symbol(fv)][fvi], mc.data[fv][fvi]; force=true)
                    end
                end
            end
        end
    end
    return nothing
end