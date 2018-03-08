module Risbo
export risbo

using ArgCheck

using WignerFunctions: Index, Cache, initial_wigner
using OffsetArrays: OffsetArray

function nosym_recurrence end

nosym(β::AbstractFloat, index::Index) = nosym(β, index, Cache{typeof(β)}())
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
    const coshb = -cos(β / 2);
    const sinhb = sin(β / 2);
    const T = promote_type(typeof(β), valtype(cache))

    intermediate = fill!(OffsetArray(T, -l:l,  -l:l), 0)

    for k ∈ -(l - 1):(l - 1), i ∈ -(l - 1):(l - 1)
        rec = nosym(β, Index(l - 1, k, i), cache) / (2l - 1)
        intermediate[i, k] += √((l - i) * (l - k)) * rec * coshb
        intermediate[i + 1, k] -= √((l + i) * (l - k)) * rec * sinhb
        intermediate[i, k + 1] += √((l - i) * (l + k)) * rec * sinhb
        intermediate[i + 1, k + 1] += √((l + i) * (l + k)) * rec * coshb
    end

    for k ∈ -l:l, i ∈ -l:l
        cache[Index(l, k, i)] = 0
    end

    for k ∈ -l:(l - 1), i ∈ -l:(l - 1)
        rec = intermediate[i + 1, k + 1] / 2l
        cache[Index(l, k, i)] += √((l - i) * (l - k)) * rec * coshb
        cache[Index(l, k, i + 1)] -= √((l + i + 1) * (l - k)) * rec * sinhb
        cache[Index(l, k + 1, i)] += √((l - i) * (l + k + 1)) * rec * sinhb
        cache[Index(l, k + 1, i + 1)] += √((l + i + 1) * (l + k + 1)) * rec * coshb
    end
end

function nosym(β::AbstractFloat)
    cache = Cache{typeof(β)}()
    (l::Int64, m₁::Int64, m₂::Int64) -> nosym(β, Index(l, m₁, m₂), cache)
end
nosym(β::Integer) = risbo(float(β))

end
