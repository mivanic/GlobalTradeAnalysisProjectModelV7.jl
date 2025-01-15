function aggComb(arr, weight, dimmap)
    uniqueNames = map(unique, dimmap)
    uniqueSizes = map(length, uniqueNames)
    if length(uniqueNames) != 1
        out = NamedArray(zeros(uniqueSizes...), (uniqueNames))
        wgt = NamedArray(zeros(uniqueSizes...), (uniqueNames))
    else
        out = NamedArray(zeros(uniqueSizes...), (uniqueNames...))
        wgt = NamedArray(zeros(uniqueSizes...), (uniqueNames...))
    end
    for i in eachindex(view(out, names(out)...))
        ii = map((j, k) -> names(out)[k][j], Tuple(i), 1:length(i))

        iii = map((j, k) -> names(dimmap[k][dimmap[k].==j]), ii, 1:length(dimmap))
        iiii = map(j -> j[1], iii)

        out[ii...] = sum(arr[iiii...] .* weight[iiii...]) / sum(weight[iiii...])
    end
    return (out)
end
