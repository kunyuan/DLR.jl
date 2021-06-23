using LinearAlgebra
using Roots
using Quadmath
# using Gaston

# const Float = Float64
const Float = BigFloat
# const Float = Float128
const Vec = Vector{Tuple{Float,Float}}
const atol = 1.0e-10
# const Λ0 = Float(1)
# const Λ = Float(100)

# """
# \int_1^∞ e^{-ω_1 τ}*e^{-ω_2*τ} dτ = 1/(ω_1+ω_2)
# """
function proj(ω1::Float, ω2::Float)
    return 1 / (ω1 + ω2)
end

function proj(freq, q, ω::Float)
    sum = Float(0)
    for (idx, ωi) in enumerate(freq)
        sum += q[idx] * proj(ωi, ω)
    end
    return sum
end

function proj(freq, q1::Vector{Float}, q2::Vector{Float})
    # println(q1)
    # println(q2)
    sum = Float(0.0)
    for (idx, ωi) in enumerate(freq)
        sum += q2[idx] * proj(freq, q1, freq[idx])
    end
    return sum
end

Norm(ω) = sqrt(proj(ω, ω))

function Norm(freq, Q, ω::Float)
#   qi=sum_j c_ij e^{-ω_j*τ}
#   norm2 = 1/2ω- \sum_i c_ij*c_ik/(ω+ω_j)/(ω+ω_k)
    norm2 = proj(ω, ω)
    # if ω > 999.0
    #     println(norm2)
    # end
    for (qi, q) in enumerate(Q)
        norm2 -= (proj(freq, q, ω))^2
        # if ω > 999.0
        #     println(norm2, " from ", proj(freq, q, ω))
        # end
    end
    # @assert norm2 > 0 "$norm2 <=0"
    norm = sqrt(abs(norm2))
end

"""
First derivative of norm2 for the vector: e^{-ωτ}-sum_i <q_i, e^{-ωτ}>
"""
function DNorm2(freq, Q, ω::Float)
#   qi=sum_j c_ij e^{-ω_j*τ}
#   dnorm2 = -1/2ω^2 + \sum_{jk} (\sum_i c_ij*c_ik)*(2ω+ω_j+ω_k)/(ω+ω_j)^2/(ω+ω_k)^2
    dnorm2 = -1 / (2 * ω^2)
    # println("omega:", ω, ", ", proj(ω, ω))
    for j in 1:length(Q)
        for k in 1:length(Q)
            amp = Float(0)
            for i in 1:length(Q)
                amp += Q[i][j] * Q[i][k]
            end
            dnorm2 += amp * (2 * ω + freq[j] + freq[k]) / (ω + freq[j])^2 / (ω + freq[k])^2
        end
    end
    return dnorm2
end

# function NormExpr(freq, Q)
#     coeff = [1.0, ]

# end

function projectedKernel(freq, Q, ω::Float)
    # K(τ, ω) ≈ \sum_i <e^{-ωτ}, qi> qi = \sum_k (\sum_ij c_ij*c_ik/(ω+ω_j) e^{-ω_k*τ})
    amp = zeros(Float, length(Q))
    for k in 1:length(Q)
        amp[k] = Float(0)
        for i in 1:length(Q)
            for j in 1:length(Q)
                amp[k] += Q[i][j] * Q[i][k] / (ω + freq[j])
            end
    end
end
    return amp
end

function orthognalize(freq, Q, ω::Float)
    idx = length(Q) + 1
    qnew = zeros(Float, length(freq))
    # println("compare ", length(freq), " vs ", length(Q[1]))
    for (qi, q) in enumerate(Q)
        qnew .-= proj(freq, q, ω) .* q
    end

    norm = Norm(freq, Q, ω)
    qnew /= norm

    return qnew, 1 / norm
end

function testOrthgonal(freq, Q)
    # maxerr = 0.0
    # pos = [1, 1]
    err = zeros(Float, (length(Q), length(Q)))
    for (i, qi) in enumerate(Q)
        for (j, qj) in enumerate(Q)
            err[i, j] = proj(freq, qi, qj)
        end
    end
    maxerr = maximum(abs.(err - I))
    println("Max Err: ", maxerr)
    @assert maxerr < atol
end

function addFreq!(freq, Q, ω)
    idx = findall(x -> x > ω, freq)
    if length(idx) == 0
        idx = length(freq) + 1
    else
        idx = idx[1]
    end
    q, q00 = orthognalize(freq, Q, ω)
    insert!(freq, idx, ω)
    insert!(Q, idx, q)
    # println("freq length = $(length(freq))")
    for qi in 1:length(Q)
        insert!(Q[qi], idx, Float(0))
        # println("Q[i] length = $(length(Q[qi]))")
    end
    Q[idx][idx] = q00
    return idx
end

function findFreqMax(freq, Q, idx, Λ)
    if idx == length(freq)
        ω = find_zero(x -> DNorm2(freq, Q, x), (freq[end] * (1 + 1e-10), freq[end] * 100), Bisection(), rtol=1e-5)
        return ω >= Λ ? Λ : ω
    else
        # println("DNorm2: ", DNorm2(freq, Q, freq[idx] * (2 + 1e-1)), " -> ", DNorm2(freq, Q, freq[idx + 1] * (0.5 - 1e-1)))
        d1, d2 = freq[idx] * (1 + 1e-5), freq[idx + 1] * (1 - 1e-5)
        if sign(DNorm2(freq, Q, d1)) == sign(DNorm2(freq, Q, d2))
            println("warning: $(freq[idx]) -> $(freq[idx + 1]) derivatives have the same sign $(DNorm2(freq, Q, d1)) -> $(DNorm2(freq, Q, d1)) !")
            ω = sqrt(freq[idx] * freq[idx + 1])
        else
            ω = find_zero(x -> DNorm2(freq, Q, x), (d1, d2), Bisection(), rtol=1e-5)
        end
        return ω
    end
end

function findFreqMedian(freq, Q, idx, Λ)
    if idx == length(freq)
        return sqrt(freq[idx] * Λ)
    else
        return sqrt(freq[idx] * freq[idx + 1])
    end
end

function scheme1(eps, Λ)
    freq = Vector{Float}([1, ])
    Q = [[1 / Norm(freq[1]), ], ]
    residual = 1.0
    candidates = [findFreqMax(freq, Q, 1, Λ), ]
    while residual > eps
        maxR = Float(0)
        idx, ifreq, newω = 1, 1, 1
        for i in 1:length(freq)
            # ω = findFreqMax(freq, Q, i)
            ω = candidates[i]
            normw = Norm(freq, Q, ω)
            if normw > maxR
                maxR = normw
                # idx = k + 1
                ifreq = i
                newω = ω 
            end
        end
        residual = maxR
        if residual > eps
            idx = addFreq!(freq, Q, newω)
            # println("add $(length(freq))")
            if idx < length(freq)
                println("$(length(freq)) basis: ω=$(Float64(newω)) between ($(Float64(freq[idx - 1])), $(Float64(freq[idx + 1])))")
            else
                println("$(length(freq)) basis: ω=$(Float64(newω)) for the last freq $(Float64(freq[idx - 1]))")
            end
            # println("residual=$residual")
            # @assert newidx == idx "idx: $idx != newidx: $newidx"
            # testOrthgonal(freq, Q)

            deleteat!(candidates, ifreq)
            @assert idx > 1 "idx=$idx"
            push!(candidates, findFreqMax(freq, Q, idx - 1, Λ))
            push!(candidates, findFreqMax(freq, Q, idx, Λ))
        end
    end
    testOrthgonal(freq, Q)
    println("residual=$residual")
return freq, Q
end

if abspath(PROGRAM_FILE) == @__FILE__    
    freq = zeros(Float, rank)
    freq[1] = 1.0
    Q = [zeros(Float, rank), ]
    Q[1][1] = 1.0 / Norm(freq[1])
    println(Q[1])
    println(Norm(freq, Q, 2.0))

    println("Derivative: ", DNorm2(freq, Q, 1.0))
    println("Derivative: ", DNorm2(freq, Q, 2.0))
    println("Derivative: ", DNorm2(freq, Q, 10.0))
    println("Derivative: ", DNorm2(freq, Q, 100.0))

    freq[2] = 2.0
    push!(Q, orthognalize(freq, Q, 2.0))
    testOrthgonal(freq, Q)

    # ω = LinRange(0.0, 20.0, 100)
    # y = [Norm(freq, Q, w) for w in ω]
    # p = plot(ω, y)
    # display(p)
    # readline()
    # println("Derivative: ", DNorm2(freq, Q, 2.0))
    # println("Derivative: ", DNorm2(freq, Q, 10.0))



end