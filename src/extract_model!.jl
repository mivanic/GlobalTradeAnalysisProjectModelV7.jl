function extract_model!(mc)
    results = merge(Dict(
            String(k) => begin
                arrayOut = NamedArray(zeros(map(length, v.axes)), v.axes)
                arrayOut[is_valid.(mc.model, v).data] .= value.(Array(v)[is_valid.(mc.model, v).data])
                arrayOut[.!is_valid.(mc.model, v).data] .= NaN
                arrayOut
            end for (k, v) in object_dictionary(mc.model)
            if v isa AbstractArray{VariableRef}
        ), Dict(
            String(k) => begin
                (is_valid(mc.model, v) ? value.(v) : NaN)
            end for (k, v) in object_dictionary(mc.model)
            if v isa VariableRef
        ))

    mc.data = merge(mc.data, Dict(k => results[k] for k ∈ setdiff(keys(results), keys(mc.parameters))))

    return nothing
end