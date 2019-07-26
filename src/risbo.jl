module Risbo

using ArgCheck

using WignerFunctions: Index, initial_wigner
using OffsetArrays: OffsetArray

const Cache = Dict{Index, T} where T <: AbstractFloat

nosym(β::AbstractFloat, index::Index) = nosym(β, index, Cache{typeof(β)}())
function nosym(β::AbstractFloat)
    cache = Cache{typeof(β)}()
    (l::Int64, m₁::Int64, m₂::Int64) -> nosym(β, Index((l, m₁, m₂)), cache)
end
nosym(β::Integer) = nosym(float(β))

function nosym(β::AbstractFloat, index::Index, cache::Cache)
    @argcheck 0 <= index.l
    @argcheck abs(index.m₁) <= index.l
    @argcheck abs(index.m₂) <= index.l

    index ∈ keys(cache) && return cache[index]

    if index.l <= 1
        result = initial_wigner(convert(valtype(cache), β), index)
        cache[index] = result
    else
        nosym_recurrence(convert(valtype(cache), β), index.l, cache)
        result = cache[index]
    end

    return result
end

function nosym_recurrence(β::AbstractFloat, l::Integer, cache::Cache)
    coshb = -cos(β / 2);
    sinhb = sin(β / 2);
    T = promote_type(typeof(β), valtype(cache))

    intermediate = fill!(OffsetArray{T}(undef, -l:l,  -l:l), 0)

    for k ∈ -(l - 1):(l - 1), i ∈ -(l - 1):(l - 1)
        rec = nosym(β, Index((l - 1, k, i)), cache) / (2l - 1)
        intermediate[i, k] += √((l - i) * (l - k)) * rec * coshb
        intermediate[i + 1, k] -= √((l + i) * (l - k)) * rec * sinhb
        intermediate[i, k + 1] += √((l - i) * (l + k)) * rec * sinhb
        intermediate[i + 1, k + 1] += √((l + i) * (l + k)) * rec * coshb
    end

    for k ∈ -l:l, i ∈ -l:l
        cache[Index((l, k, i))] = 0
    end

    for k ∈ -l:(l - 1), i ∈ -l:(l - 1)
        rec = intermediate[i + 1, k + 1] / 2l
        cache[Index((l, k, i))] += √((l - i) * (l - k)) * rec * coshb
        cache[Index((l, k, i + 1))] -= √((l + i + 1) * (l - k)) * rec * sinhb
        cache[Index((l, k + 1, i))] += √((l - i) * (l + k + 1)) * rec * sinhb
        cache[Index((l, k + 1, i + 1))] += √((l + i + 1) * (l + k + 1)) * rec * coshb
    end
end


half(β::AbstractFloat, index::Index) = half(β, index, Cache{typeof(β)}())
function half(β::AbstractFloat)
    cache = Cache{typeof(β)}()
    (l::Int64, m₁::Int64, m₂::Int64) -> half(β, Index((l, m₁, m₂)), cache)
end
half(β::Integer) = half(float(β))

function half(β::AbstractFloat, index::Index, cache::Cache)
    @argcheck 0 <= index.l
    @argcheck abs(index.m₁) <= index.l
    @argcheck abs(index.m₂) <= index.l

    if index.m₂ <= 0
        index ∈ keys(cache) && return cache[index]

        if index.l <= 1
            cache[index] = initial_wigner(convert(valtype(cache), β), index)
        else
            half_recurrence(convert(valtype(cache), β), index.l, cache)
            cache[index]
        end
    elseif (index.m₁ + index.m₂) % 2 == 0
        half(β, Index((index.l, -index.m₁, -index.m₂)), cache)
    else
        -half(β, Index((index.l, -index.m₁, -index.m₂)), cache)
    end
end

function half_recurrence(β::AbstractFloat, l::Integer, cache::Cache)
    coshb = -cos(β / 2);
    sinhb = sin(β / 2);
    T = promote_type(typeof(β), valtype(cache))

    intermediate = fill!(OffsetArray(T, -l:l,  -l:l), 0)

    for k ∈ -(l - 1):(l - 1), i ∈ -(l - 1):1
        rec = half(β, Index((l - 1, k, i)), cache) / (2l - 1)
        intermediate[i, k] += √((l - i) * (l - k)) * rec * coshb
        intermediate[i + 1, k] -= √((l + i) * (l - k)) * rec * sinhb
        intermediate[i, k + 1] += √((l - i) * (l + k)) * rec * sinhb
        intermediate[i + 1, k + 1] += √((l + i) * (l + k)) * rec * coshb
    end

    for k ∈ -l:l, i ∈ -l:0
        cache[Index((l, k, i))] = 0
    end

    for k ∈ -l:(l - 1), i ∈ -l:0
        rec = intermediate[i + 1, k + 1] / 2l
        cache[Index((l, k, i))] += √((l - i) * (l - k)) * rec * coshb
        cache[Index((l, k, i + 1))] -= √((l + i + 1) * (l - k)) * rec * sinhb
        cache[Index((l, k + 1, i))] += √((l - i) * (l + k + 1)) * rec * sinhb
        cache[Index((l, k + 1, i + 1))] += √((l + i + 1) * (l + k + 1)) * rec * coshb
    end
end

end
