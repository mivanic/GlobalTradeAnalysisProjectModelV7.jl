function get_sample_data()
    return JLD2.load_object(joinpath(@__DIR__, ".", "assets", "example.jld2"))
end