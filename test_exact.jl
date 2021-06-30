include("./src/exact.jl")

using Test

@testset "Correlator Representation" begin
    rtol(x, y) = maximum(abs.(x - y)) / maximum(abs.(x))
    Λ = 10000
    eps = 1.0e-8
    freq, Q = scheme1(eps, Λ)

    t = 

end