export histogram

"""
    histogram(array; bins=minimum([100, length(array)]), normalization=:uniform, custom_range=false)

# Description
    Compute the histogram of an array. Useful for barplot in GLMakie.

# Arguments
- `array::AbstractArray`: Array to compute the histogram of.

# Keyword Arguments
- `bins::Integer`: Number of bins to use.
- `normalization::AbstractArray`: Normalization to use. If `:uniform`, then the normalization is uniform.
- `custom_range::Tuple`: Custom range to use. If `false`, then the range is computed from the data.

# Returns
- `Tuple{AbstractArray, AbstractArray}`: Tuple of the bin centers and the histogram values.
"""
function histogram(
    array;
    bins=minimum([100, length(array)]),
    normalization=:uniform,
    custom_range=false
)
    tmp = zeros(bins)
    if custom_range isa Tuple
        down, up = custom_range
    else
        down, up = extrema(array)
    end
    down, up = down == up ? (down - 1, up + 1) : (down, up) # edge case
    bucket = collect(range(down, up, length=bins + 1))
    if normalization == :uniform
        normalization = ones(length(array)) ./ length(array)
    end
    for i in eachindex(array)
        # normalize then multiply by bins
        val = (array[i] - down) / (up - down) * bins
        ind = ceil(Int, val)
        # handle edge cases
        ind = maximum([ind, 1])
        ind = minimum([ind, bins])
        tmp[ind] += normalization[i]
    end
    return (bucket[2:end] + bucket[1:(end-1)]) .* 0.5, tmp
end
