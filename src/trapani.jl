module Trapani

using ArgCheck
using MicroLogging
using NamedTuples: @NT

using WignerFunctions: Index, initial_wigner
using OffsetArrays: OffsetArray

const CacheIndex = @NT(branch::Bool, l::Int64, m₁::Int64, m₂::Int64)
const Cache = Dict{CacheIndex, T} where T <: AbstractFloat
function Base.getindex(cache::Cache, branch::Bool, l::Integer, m₁::Integer, m₂::Integer)
    cache[CacheIndex(branch, l, m₁, m₂)]
end
function Base.setindex!(cache::Cache, value::Number, branch::Bool, l::Integer, m₁::Integer, m₂::Integer)
    cache[CacheIndex(branch, l, m₁, m₂)] = value
end
function Base.getindex(cache::Cache, branch::Bool, index::Index)
    cache[CacheIndex(branch, index.l, index.m₁, index.m₂)]
end
function Base.setindex!(cache::Cache, value::Number, branch::Bool, index::Index)
    cache[CacheIndex(branch, index.l, index.m₁, index.m₂)] = value
end


trapani(β::AbstractFloat, index::Index) = trapani(β, index, Cache{typeof(β)}())
function trapani(β::AbstractFloat)
    cache = Cache{typeof(β)}()
    (l::Int64, m₁::Int64, m₂::Int64) -> trapani(β, Index(l, m₁, m₂), true, cache)
end
trapani(β::Integer) = trapani(float(β))
function trapani(β::AbstractFloat, l::Integer, m₁::Integer, m₂::Integer, branch::Bool, cache::Cache)
    trapani(β, Index(l, m₁, m₂), branch, cache)
end

function trapani(β::AbstractFloat, index::Index, branch, cache::Cache)
    @argcheck 0 <= index.l
    @argcheck abs(index.m₁) <= index.l
    @argcheck abs(index.m₂) <= index.l

    factor(n::Integer) = x -> n % 2 == 0 ? x: -x
    if index.m₂ < 0 
        @debug "m₂ < 0" branch β index
        trapani(β, index.l, index.m₁, -index.m₂, !branch, cache) |> factor(index.l + index.m₁)
    elseif index.m₁ < 0
        @debug "m₁ < 0" branch β index
        trapani(β, index.l, -index.m₁, index.m₂, !branch, cache) |> factor(index.l + index.m₂)
    elseif index.m₂ > index.m₁
        @debug "m₂ > m₁" branch β index
        trapani(β, index.l, index.m₂, index.m₁, branch, cache) |> factor(index.m₁ + index.m₂)
    elseif CacheIndex(branch, index.l, index.m₁, index.m₂) ∈ keys(cache)
        @debug "cache" branch β index
        cache[branch, index]
    elseif index.l <= 1
        if branch
            cache[branch, index] = initial_wigner(convert(valtype(cache), β), index)
        else
            cache[branch, index] = initial_wigner(convert(valtype(cache), β - π), index)
        end
        @debug "init" branch β index value=cache[branch, index]
        cache[branch, index]
    else
        recurrence(convert(valtype(cache), β), index.l, branch, cache)
        @debug "rec end" branch β index value=cache[branch, index]
        cache[branch, index]
    end
end

function recurrence(β::AbstractFloat, l::Integer, branch::Bool, cache::Cache)
    recursion = (l, m₁, m₂) -> trapani(β, l, m₁, m₂, branch, cache)

    cache[branch, l, l, 0] = -√((2l - 1) / 2l) * recursion(l - 1, l - 1, 0)
    for m₂ ∈ 1:l
        factor = √(l * (2l - 1) / (2(l + m₂) * (l + m₂ - 1)))
        cache[branch, l, l, m₂] = factor * recursion(l - 1, l - 1, m₂ - 1)
        @debug "rec: " l m₂ factor cache[branch, l, l, m₂]
    end

    t1 = (l, m₁, m₂) -> 2m₂ / √((l - m₁) * (l + m₁ + 1))
    t2 = (l, m₁, m₂) -> √((l - m₁ - 1) * (l + m₁ + 2) / (l - m₁) * (l + m₁ + 1))
    for m₂ ∈ 0:l
        cache[branch, l, l - 1, m₂] = 2m₂ / √(2l) * recursion(l, l, m₂ - 1)
        for m₁ ∈ (l - 2):-1:m₂
            a = t1(l, m₁, m₂) * recursion(l, m₁ + 1, m₂)
            b = t2(l, m₁, m₂) * recursion(l, m₁ + 2, m₂)
            cache[branch, l, m₁, m₂] = a - b
        end
    end
end

end
