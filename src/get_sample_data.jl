"""
    get_sample_data()

Returns a small sample data set from publicly available version 9 of the GTAP database.

### Args:
- None

### Retruns:
A named tuple with
- `hData`: a dictionary with sample aggregated data
- `hParamters`: a dictionary with sample aggregated parameters
- `hSets`: a dictionary with sample aggregated sets
       
### Example:

`julia
(; hData, hParameters, hSets) = get_sample_data()
`
"""
function get_sample_data()
    return JLD2.load_object(joinpath(@__DIR__, ".", "assets", "example.jld2"))
end