module WignerFunctions
export risbo

using NamedTuples: @NT
using ArgCheck

const Index = @NT(l::Int64, m₁::Int64, m₂::Int64)
const Cache = Dict{Index, T} where T <: AbstractFloat

function initial_wigner(β::Number, index::Index)
    @argcheck 0 <= index.l < 2
    @argcheck abs(index.m₁) <= index.l
    @argcheck abs(index.m₂) <= index.l

    if index.l == 0
        typeof(β)(1)
    elseif index.m₁ == -1 && index.m₂ == -1
        cos(β/2) * cos(β/2)
    elseif index.m₁ == -1 && index.m₂ == 0
        sin(β) / √(typeof(β)(2))
    elseif index.m₁ == -1 && index.m₂ == 1
        sin(β/2) * sin(β/2)
    elseif index.m₁ == 0 && index.m₂ == -1
        -sin(β) / √(typeof(β)(2))
    elseif index.m₁ == 0 && index.m₂ == 0
        cos(β)
    elseif index.m₁ == 0 && index.m₂ == 1
        sin(β) / √(typeof(β)(2))
    elseif index.m₁ == 1 && index.m₂ == -1
        sin(β/2) * sin(β/2)
    elseif index.m₁ == 1 && index.m₂ == 0
        -sin(β) / √(typeof(β)(2))
    elseif index.m₁ == 1 && index.m₂ == 1
        cos(β/2) * cos(β/2)
    end
end

include("risbo.jl")
risbo = Risbo.nosym

function naive(β::AbstractFloat, index::Index)
    result::typeof(β) = 0
    j::BigInt, m::BigInt, n::BigInt = index
    numerator = √(factorial(j + m) * factorial(j - m) * factorial(j + n) * factorial(j - n))
    for k ∈ max(0, n - m):min(j - m, j + n)
        denumerator = (
            factorial(j - m - k) * factorial(j + n - k) * factorial(k + m - n)
            * factorial(k)
        )
        factor = ((k + m - n) % 2 == 0 ? 1: -1) * numerator / denumerator
        result += factor * cos(β/2)^(2j - 2k + n - m) * sin(β/2)^(2k + m - n)
    end
    result
end
naive(β::AbstractFloat, l::Integer, m₁::Integer, m₂::Integer) = naive(β, Index(l, m₁, m₂))

end # module