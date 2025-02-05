"""
    extract_model!(model_container)
"""
function extract_model!(model_container)
    results = merge(Dict(
            String(k) => begin
                arrayOut = NamedArray(zeros(map(length, v.axes)), v.axes)
                arrayOut[is_valid.(model_container.model, v).data] .= value.(Array(v)[is_valid.(model_container.model, v).data])
                arrayOut[.!is_valid.(model_container.model, v).data] .= NaN
                arrayOut
            end for (k, v) in object_dictionary(model_container.model)
            if v isa AbstractArray{VariableRef}
        ), Dict(
            String(k) => begin
                (is_valid(model_container.model, v) ? value.(v) : NaN)
            end for (k, v) in object_dictionary(model_container.model)
            if v isa VariableRef
        ))

    model_container.data = merge(model_container.data, Dict(k => results[k] for k ∈ setdiff(keys(results), keys(model_container.parameters))))
    return nothing
end