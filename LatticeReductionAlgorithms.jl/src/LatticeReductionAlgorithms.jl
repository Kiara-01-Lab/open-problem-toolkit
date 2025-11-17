module LatticeReductionAlgorithms

using LinearAlgebra
using OffsetArrays

using LinearAlgebra
using OffsetArrays: OffsetVector

struct GSOData{I<:Integer,F<:Real}
    B::Matrix{I}
    B⃗::Vector{F}
    Q::Matrix{F}
    R::Matrix{F}
end

function GSOData(B::AbstractMatrix)
    F = qr(B)
    R̃ = F.R
    for i in axes(R̃, 1)
        dᵢ = R̃[i, i]
        for j = i:size(R̃, 2)
            R̃[i, j] /= dᵢ
        end
    end
    Q̃ = F.Q * Diagonal(F.R)
    B⃗ = [dot((@view Q̃[:, j]), (@view Q̃[:, j])) for j in axes(Q̃, 2)]
    return GSOData(Matrix(B), B⃗, Q̃, R̃)
end

function partial_size_reduce!(g::GSOData{IType,FType}, i::Int, j::Int) where {IType,FType}
    if !(i < j)
        throw(ArgumentError("should satisfy i < j, actual i=$(i), j=$(j)"))
    end
    μ_ij = g.R[i, j]
    q = fplll_round(IType, μ_ij)
    negq = -q
    bi = @view g.B[:, i]
    bj = @view g.B[:, j]
    @inbounds for r in eachindex(bi, bj)
        bj[r] = muladd(negq, bi[r], bj[r])
    end
    @inbounds for l = 1:i
        g.R[l, j] = muladd(negq, g.R[l, i], g.R[l, j])
    end
    return g
end

@inline function iszerovec(v)
    r = true
    for i in eachindex(v)
        r *= v[i] ≈ zero(v[i])
    end
    r
end

function fplll_round(T, μ)
    # fplll の最近整数関数に合わせる: 最近接（0.5 はゼロから遠い側に丸め）
    if μ >= zero(μ)
        q_i = floor(T, μ + (one(μ) / 2))
    else
        q_i = ceil(T, μ - (one(μ) / 2))
    end
end


function partial_size_reduce!(
    B::AbstractMatrix{T},
    R::AbstractMatrix,
    i::Int,
    k::Int,
) where {T}
    (1 ≤ i < k) || return
    μ = R[i, k]
    m = fplll_round(T, μ)
    m == 0 && return
    B[:, k] .-= m .* B[:, i]
    @inbounds begin
        for j = 1:(i-1)
            R[j, k] -= m * R[j, i]
        end
        R[i, k] -= m
    end
end

function size_reduce!(g::GSOData)
    R = g.R
    for j = 2:size(R, 2)
        for i = (j-1):-1:1
            partial_size_reduce!(g, i, j)
        end
    end
    g
end

function ENUM_reduce(
    μ::AbstractMatrix{T},
    B⃗::AbstractVector{T},
    R²::AbstractVector,
)::Tuple{Vector{BigInt},Bool} where {T<:AbstractFloat}
    n = length(B⃗)
    σ = zeros(T, n+1, n)
    r = OffsetVector(collect(0:n), 0:n)
    ρ = zeros(T, n+1)
    v = zeros(BigInt, n)
    v[begin] = 1
    c = zeros(T, n)
    w = zeros(BigInt, n)
    last_nonzero = 1
    k = 1
    while true
        ρ[k] = ρ[k+1] + (v[k] - c[k]) ^ 2 * B⃗[k]
        if ρ[k] ≤ R²[n+1-k]
            if k == 1
                # return solution
                return (v, true)
            end
            k -= 1
            r[k-1] = max(r[k-1], r[k])
            for i = r[k]:-1:(k+1)
                σ[i, k] = σ[i+1, k] + μ[k, i] * v[i]
            end
            c[k] = -σ[k+1, k]
            v[k] = fplll_round(T, c[k])
            w[k] = 1
        else
            k += 1
            if k == n + 1
                # solution not found
                return (zeros(n), false)
            end
            r[k-1] = k
            if k ≥ last_nonzero
                last_nonzero = k
                v[k] += 1
            else
                if v[k] > c[k]
                    v[k] -= w[k]
                else
                    v[k] += w[k]
                end
                w[k] += 1
            end # if
        end # if
    end # while
end # function

function gsoupdate!(g::GSOData, k::Integer)
    k ≥ 2 || throw(ArgumentError("k should satisfy k ≥ 2, actual k=$(k)"))
    for i in axes(g.B, 1)
        # swap
        g.B[i, k-1], g.B[i, k] = g.B[i, k], g.B[i, k-1]
    end
    μ = g.R
    ν = μ[k-1, k]
    B = g.B⃗[k] + ν ^ 2 * g.B⃗[k-1]
    μ[k-1, k] = ν * g.B⃗[k-1] / B
    g.B⃗[k] = g.B⃗[k] * g.B⃗[k-1] / B
    g.B⃗[k-1] = B
    for j = 1:(k-2)
        # swap
        μ[j, k-1], μ[j, k] = μ[j, k], μ[j, k-1]
    end
    n = size(μ, 2)
    for i = (k+1):n
        t = μ[k, i]
        μ[k, i] = μ[k-1, i] - ν * t
        μ[k-1, i] = t + μ[k-1, k] * μ[k, i]
    end
    g
end

function LLL_reduce(B::AbstractMatrix, δ::Real)
    if !(0.25 < δ < 1)
        throw(ArgumentError("Input δ must satisfy 0.25 < δ < 1"))
    end
    g = GSOData(B)
    k = 2
    n = size(g.B, 2)
    while k ≤ n
        for j = (k-1):-1:1
            partial_size_reduce!(g, j, k)
        end
        if g.B⃗[k] ≥ (δ - g.R[k-1, k] ^ 2) * g.B⃗[k-1]
            # Lovász 条件を満たす
            k += 1
        else
            # swap basis
            gsoupdate!(g, k)
            k = max(k-1, 2)
        end
    end
    return g
end

function _MLLL_reduce!(ℬ::AbstractMatrix{T}, δ::Float64) where {T}
    h = size(ℬ, 2)
    z = h
    g = 1

    # GSO 用のワーク
    Q = float(T).(ℬ)
    R = zeros(eltype(Q), size(ℬ))
    for i = 1:min(size(R, 1), size(R, 2))
        R[i, i] = 1
    end
    B⃗ = zeros(eltype(Q), size(ℬ, 2))

    while g ≤ z
        b_g = @view ℬ[:, g]

        # 0 列なら末尾と交換して z を詰める
        if iszerovec(b_g)
            if g < z
                v = ℬ[:, g]
                ℬ[:, g] .= ℬ[:, z]
                ℬ[:, z] .= v
            end
            z -= 1
        end

        # GSO: b_g* を直交化
        Q[:, g] .= ℬ[:, g]
        for i = 1:(g-1)
            b_i_ast = @view Q[:, i]
            B_i = dot(b_i_ast, b_i_ast)
            μ_ig = dot(Q[:, g], b_i_ast) / B_i
            R[i, g] = μ_ig
            Q[:, g] .-= μ_ig .* b_i_ast
        end
        B⃗[g] = dot((@view Q[:, g]), (@view Q[:, g]))
        if g == 1
            g = 2
            continue
        end

        # --- MLLL 本体 ---
        l = g
        k = g
        startagain = false
        while (k ≤ l) && !startagain
            partial_size_reduce!(ℬ, R, k - 1, k)

            ν = R[k-1, k]
            B = B⃗[k] + ν^2 * B⃗[k-1]

            if B ≥ δ * B⃗[k-1]
                for j = (k-2):-1:1
                    partial_size_reduce!(ℬ, R, j, k)
                end
                k += 1
            else
                if iszerovec(@view ℬ[:, k])
                    if k < z
                        v = ℬ[:, k]
                        ℬ[:, k] .= ℬ[:, z]
                        ℬ[:, z] .= v
                    end
                    z -= 1
                    g = k
                    startagain = true
                else
                    ℬ[:, k-1], ℬ[:, k] = ℬ[:, k], ℬ[:, k-1]
                    for j = 1:(k-2)
                        R[j, k], R[j, k-1] = R[j, k-1], R[j, k]
                    end

                    if !(B ≈ 0)
                        if B⃗[k] ≈ 0
                            B⃗[k] = B
                            Q[:, k-1] .= ν .* Q[:, k-1]
                            R[k-1, k] = inv(ν)
                            for i = (k+1):l
                                R[k-1, i] /= ν
                            end
                        else
                            t = B⃗[k-1] / B
                            R[k-1, k] = ν * t
                            w = Q[:, k-1]
                            Q[:, k-1] .= Q[:, k] .+ ν .* w
                            B⃗[k-1] = B
                            if k ≤ l
                                Q[:, k] .= -R[k-1, k] .* Q[:, k] .+ (B⃗[k] / B) .* w
                                B⃗[k] *= t
                            end
                            for i = (k+1):l
                                t = R[k, i]
                                R[k, i] = R[k-1, i] - ν * t
                                R[k-1, i] = t + R[k-1, k] * R[k, i]
                            end
                        end
                    else
                        B⃗[k], B⃗[k-1] = B⃗[k-1], B⃗[k]
                        Q[:, k], Q[:, k-1] = Q[:, k-1], Q[:, k]
                        for i = (k+1):l
                            R[k, i], R[k-1, i] = R[k-1, i], R[k, i]
                        end
                    end
                    k = max(k - 1, 2)
                end
            end
        end

        if !startagain
            g += 1
        end
    end
end

function MLLL_reduce!(B, δ)
    g_target_j = typemax(Int)
    while true
        _MLLL_reduce!(B, δ)
        target_j = -1
        for j = size(B, 2):-1:1
            if iszerovec(@view B[:, j])
                target_j = j
            end
        end

        if target_j < 0
            break
        end
        if target_j == g_target_j
            break
        end

        _MLLL_reduce!((@view B[:, 1:target_j]), δ)
        g_target_j = target_j
        break
    end
    B
end

function find_svp_by_enum(B, k, l)
    @assert k < l
    ε = 0.99
    g = GSOData(B)
    n = size(B, 2)

    R²ₙ = ε * maximum(g.B⃗[k:l])
    R² = [R²ₙ for _ = k:l]

    μ = g.R[k:l, k:l]
    B⃗ = g.B⃗[k:l]
    coeff, is_succeeded = ENUM_reduce(μ, B⃗, R²)
    return coeff
end

function BKZ_reduce!(B::AbstractMatrix, β::Integer, δ::Real)
    g = LLL_reduce(B, δ)
    B .= g.B
    n = size(B, 2)
    z = 0
    k = 0
    while z < n - 1
        k = mod(k, n-1) + 1
        l = min(k + β - 1, n)
        h = min(l + 1, n)
        coeff = find_svp_by_enum(B, k, l)
        if iszerovec(coeff)
            # @info "coeff got zero vector"
            break
        end
        v = zeros(eltype(B), n)
        for (i, idx_k_to_l) in enumerate(k:l)
            v += coeff[i] * B[:, idx_k_to_l]
        end

        g = GSOData(B)

        πₖv_norm2 = 0
        for i = k:n
            πₖv_norm2 += dot(v, @view(g.Q[:, i])) ^ 2 / dot(@view(g.Q[:, i]), @view(g.Q[:, i]))
        end
        if norm(g.Q[:, k]) > sqrt(πₖv_norm2) + 0.00001
            z = 0
            # Bsub = hcat(((@view B[:, i]) for i = 1:(k-1))..., v, ((@view B[:, i]) for i = k:h)...)
            Bsub = similar(B, size(B, 1), h+1)
            for i = 1:(k-1)
                Bsub[:, i] .= @view(B[:, i])
            end
            Bsub[:, k] .= v
            for (i, idx) in enumerate(k:h)
                Bsub[:, k+i] .= @view(B[:, idx])
            end
            MLLL_reduce!(Bsub, δ)
            B[:, 1:h] .= @view Bsub[:, 1:h]
        else
            z += 1
            g′ = LLL_reduce((@view B[:, 1:h]), δ)
            B[:, 1:h] = g′.B
        end
    end # while
end
end # module LatticeReductionAlgorithms
