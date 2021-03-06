using WignerFunctions
using Test
using Random


d222(β::AbstractFloat) = (1 + cos(β))^2 / 4
d221(β::AbstractFloat) = -sin(β) * (1 + cos(β)) / 2
d220(β::AbstractFloat) = sin(β) * sin(β) * √(3 / 8)
d221̄(β::AbstractFloat) = -sin(β) * (1 - cos(β)) / 2
d222̄(β::AbstractFloat) = (1 - cos(β))^2 / 4
d211(β::AbstractFloat) = (2cos(β) * cos(β) + cos(β) - 1) / 2
d210(β::AbstractFloat) = -sin(2β) * √(3 / 8)
d211̄(β::AbstractFloat) = (-2cos(β) * cos(β) + cos(β) + 1) / 2
d200(β::AbstractFloat) = (3cos(β) * cos(β) - 1) / 2
WignerFunctions.Index(::typeof(d222)) = WignerFunctions.Index((2, 2, 2))
WignerFunctions.Index(::typeof(d221)) = WignerFunctions.Index((2, 2, 1))
WignerFunctions.Index(::typeof(d220)) = WignerFunctions.Index((2, 2, 0))
WignerFunctions.Index(::typeof(d221̄)) = WignerFunctions.Index((2, 2, -1))
WignerFunctions.Index(::typeof(d222̄)) = WignerFunctions.Index((2, 2, -2))
WignerFunctions.Index(::typeof(d211)) = WignerFunctions.Index((2, 1, 1))
WignerFunctions.Index(::typeof(d210)) = WignerFunctions.Index((2, 1, 0))
WignerFunctions.Index(::typeof(d211̄)) = WignerFunctions.Index((2, 1, -1))
WignerFunctions.Index(::typeof(d200)) = WignerFunctions.Index((2, 0, 0))


@testset "Naive" begin
    naive = WignerFunctions.naive
    @testset "diagonal: β=-0, l=$l, m₁=$m₁, m₂=$m₂" for l ∈ 0:4, m₁ ∈ -l:l, m₂ ∈ -l:l
        @test naive(BigFloat(0), l, m₁, m₂) ≈ (m₁ == m₂ ? 1 : 0)
    end

    functions = [d222, d221, d220, d221̄, d222̄, d211, d210, d211̄, d211̄, d200]
    @testset "function $func" for func ∈ functions
        β = 2π*rand(10)
        l, m, m′ = WignerFunctions.Index(func)
        @test all(func.(β) .≈ naive.(BigFloat.(β), l, m, m′))
        @test all(func.(β) .≈ ((m - m′) % 2 == 0 ? 1 : -1) .* naive.(BigFloat.(β),  l, m′, m))
        @test all(naive.(BigFloat.(β), l, m′, m) .≈ naive.(BigFloat.(β),  l, -m, -m′))
    end
end

@testset "$method" for method ∈ (WignerFunctions.Risbo.nosym,)
    dmatrix = method(0)
    @testset "diagonal: β=-0, l=$l, m₁=$m₁, m₂=$m₂" for l ∈ 0:4, m₁ ∈ -l:l, m₂ ∈ -l:l
        @test dmatrix(l, m₁, m₂) ≈ (m₁ == m₂ ? 1 : 0)
    end

    l = rand(2:20, 50)
    m = rand.(range.(-l, l, step=1))
    m′ = rand.(range.(-l, l, step=1))
    β = 2π * rand(length(l))
    @testset "Symmetries" begin
        @testset "β=$(β[i]), l=$(l[i]), m=$(m[i]), m′=$(m′[i]) " for i ∈ 1:length(l)
            k, m₁, m₂ = l[i], m[i], m′[i]
            dmatrix = method(β[i])
            @test dmatrix(k, m₁, m₂) ≈ ((m₁ - m₂) % 2 == 0 ? 1 : -1) * dmatrix(k, m₂, m₁)
            @test dmatrix(k, m₁, m₂) ≈ ((m₁ - m₂) % 2 == 0 ? 1 : -1) * dmatrix(k, -m₁, -m₂)
            @test dmatrix(k, m₁, m₂) ≈ dmatrix(k, -m₂, -m₁)
        end
    end
 
    naive = WignerFunctions.naive
    oldprec = precision(BigFloat)
    setprecision(512)
    @testset "Against naive β=$(β[i]), l=$(l[i]), m=$(m[i]), m′=$(m′[i]) " for i ∈ 1:length(l)
        @test method(β[i])(l[i], m[i], m′[i]) ≈ float(naive(BigFloat(β[i]), l[i], m[i], m′[i]))
    end
    setprecision(oldprec)
end
nothing

@testset "Delta recurrence" begin
    recurrence = WignerFunctions.Trapani.delta_recurrence
    @test recurrence(0, 0, 0) == 1
    @test recurrence(3, 4, 0) == 0
    @test recurrence(3, 0, -4) == 0


    ls = rand(2:20, 50)
    ms = rand.(range.(-ls, ls, step=1))
    ms′ = rand.(range.(-ls, ls, step=1))
    oldprec = precision(BigFloat)
    setprecision(512)
    naive = (l, m₁, m₂) -> float(WignerFunctions.naive(BigFloat(pi)/2, l, m₁, m₂))
    @testset "Against naive l=$(ls[i]), m=$(ms[i]), m′=$(ms′[i]) " for i ∈ 1:length(ls)
        l, m₁, m₂ = ls[i], ms[i], ms′[i]
        @test recurrence(l, m₁, m₂) ≈ naive(l, m₁, m₂)
    end
    setprecision(oldprec)
end
