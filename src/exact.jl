using LinearAlgebra

const Float = Float64
const Vec = Vector{Tuple{Float,Float}}
const rank = 16
const atol = 1.0e-10
const Λ = Float(1000)

# """
# \int_1^∞ e^{-ω_1 τ}*e^{-ω_2*τ} dτ = 1/(ω_1+ω_2)
# """
function proj(ω1::Float, ω2::Float)
    return 1 / (ω1 + ω2)
end

function proj(freq, q, ω::Float)
    sum = Float(0)
    for (idx, ωi) in enumerate(freq)
        # println("idx: $idx, ωi: $ωi")
        if ωi >= 1 - eps(Float(0))
            sum += q[idx] * proj(ωi, ω)
            # println("more: $idx, ", q[idx], ", ", ωi, ", ", proj(ωi, ω))
        end
    end
    return sum
end

function proj(freq, q1::Vector{Float}, q2::Vector{Float})
    # println(q1)
    # println(q2)
    sum = Float(0.0)
    for (idx, ωi) in enumerate(freq)
        if ωi >= 1 - eps(Float(0))
            sum += q2[idx] * proj(freq, q1, freq[idx])
            # println(q2[idx], ", ", freq[idx], ", ", proj(freq, q1, freq[idx]))
        end
    end
    return sum
end

Norm(ω) = sqrt(proj(ω, ω))

function Norm(freq, Q, ω::Float)
#   qi=sum_j c_ij e^{-ω_j*τ}
#   norm2 = 1/2ω- \sum_i c_ij*c_ik/(ω+ω_j)/(ω+ω_k)
    norm2 = proj(ω, ω)
    for (qi, q) in enumerate(Q)
        norm2 -= (proj(freq, q, ω))^2
    end
    @assert norm2 > 0
    norm = sqrt(norm2)
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

function orthognalize(freq, Q, ω::Float)
    idx = length(Q) + 1
    qnew = zeros(Float, rank)
    # println(idx)
    qnew[idx] = 1.0
    # println(qnew)
    for (qi, q) in enumerate(Q)
        qnew .-= proj(freq, q, ω) .* q
    end

    norm = Norm(freq, Q, ω)
    qnew /= norm

    return qnew
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


if abspath(PROGRAM_FILE) == @__FILE__    
    freq = zeros(Float, rank)
    freq[1] = 1.0
    Q = [zeros(Float, rank), ]
    Q[1][1] = 1.0 / Norm(freq[1])
    println(Q[1])
    println(Norm(freq, Q, 2.0))

    # println("Derivative: ", DNorm2(freq, Q, 1.0))
    # println("Derivative: ", DNorm2(freq, Q, 2.0))
    # println("Derivative: ", DNorm2(freq, Q, 10.0))
    # println("Derivative: ", DNorm2(freq, Q, 100.0))

    freq[2] = 2.0
    push!(Q, orthognalize(freq, Q, 2.0))
    testOrthgonal(freq, Q)

    # println("Derivative: ", DNorm2(freq, Q, 2.0))
    # println("Derivative: ", DNorm2(freq, Q, 10.0))



end