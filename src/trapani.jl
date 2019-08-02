module Trapani

using ArgCheck

using WignerFunctions.Risbo: initial_wigner, Cache

function delta_recurrence(l::Integer, m₁::Integer, m₂::Integer, cache::Cache)
    if l < 0 || abs(l) < m₁ || abs(l) < m₂
        return 0
    elseif m₁ < 0 && m₂ < 0
        return delta_recurrence(l, -m₂, -m₁, cache)
    elseif m₁ < 0
        return ((l - m₂) % 2 == 0 ? 1 : -1) * delta_recurrence(l, -m₁, m₂, cache)
    elseif m₂ < 0
        return ((l - m₁) % 2 == 0 ? 1 : -1) * delta_recurrence(l, m₁, -m₂, cache)
    elseif m₂ < m₁
        return ((m₁ - m₂) % 2 == 0 ? 1 : -1) * delta_recurrence(l, m₂, m₁, cache)
    elseif (l=l, m₁=m₁, m₂=m₂) ∈ keys(cache)
        return cache[(l=l, m₁=m₁, m₂=m₂)]
    end
    result = if l == 0 && m₁ == 0 && m₂ == 0 
            1
        elseif l == m₁
            (
                sqrt((2l*l - l) / (2 * (l + m₂) * (l + m₂ - 1)))
                * delta_recurrence(l - 1, l - 1, m₂ - 1, cache)
            )
        else
            (
                2m₂ / sqrt((l - m₁) * (l + m₁ + 1))
                * delta_recurrence(l, m₁ + 1, m₂, cache)
                - sqrt((l - m₁ - 1) * (l + m₁ + 2) / ((l - m₁) * (l + m₁ + 1)))
                * delta_recurrence(l, m₁ + 2, m₂, cache)
            )
    end
    cache[(l=l, m₁=m₁, m₂=m₂)] = result
    result
end

function delta_recurrence(l::Integer, m₁::Integer, m₂::Integer)
    delta_recurrence(l, m₁, m₂, Cache{Float64}())
end

end
