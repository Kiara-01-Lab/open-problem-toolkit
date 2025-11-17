### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 67320cd6-a8a2-11f0-81d5-ff998c22343b
begin
    using Test
    using OffsetArrays
    using LinearAlgebra
end

# ╔═╡ bb3554a0-27c4-4579-a0c6-863a23934709
md"""
ENUM アルゴリズムを利用した格子上の最短ベクトルの数え上げ
"""

# ╔═╡ 851d1600-74b6-4049-8a8d-94638b688e93
begin
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
        return GSOData(B, B⃗, Q̃, R̃)
    end

    function partial_size_reduce!(
        g::GSOData{IType,FType},
        i::Int,
        j::Int,
    ) where {IType,FType}
        if !(i < j)
            throw(ArgumentError("should satisfy i < j, actual i=$(i), j=$(j)"))
        end
        μ_ij = g.R[i, j]
        q = round(IType, μ_ij)
        bi = @view g.B[:, i]
        bj = @view g.B[:, j]
        @. bj -= q * bi
        for l = 1:i
            g.R[l, j] -= q * g.R[l, i]
        end
        g
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

    @testset "GSO vector is invariant" begin
        B = [
            5 2 3
            -3 -7 -10
            -7 -7 0
        ]
        Q = GSOData(B).Q
        B_reduced = [
            5 -3 1
            -3 -4 -3
            -7 0 7
        ]
        Q_reduced = GSOData(B_reduced).Q
        @test isapprox(Q, Q_reduced, atol = 1e-14)
    end
end

# ╔═╡ 17033c43-3f58-4cae-9af1-83d87249019d
function enum_algorithm(
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
            v[k] = round(BigInt, c[k])
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

# ╔═╡ 0578cc9c-6620-41b3-b9e3-5304b04eea6b
function find_svp_by_enum(B)
    ε = 0.99
    g = GSOData(B)
    n = size(B, 2)
    R²ₙ = ε * norm(g.B⃗[1])
    R² = [k * R²ₙ / n for k = 1:n]

    μ = g.R
    B⃗ = g.B⃗
    v = zeros(eltype(B), n)
    while true
        coeff, is_succeeded = enum_algorithm(μ, B⃗, R²)
        if is_succeeded
            fill!(v, zero(eltype(B)))
            for i in eachindex(coeff)
                v += coeff[i] * B[:, i]
            end
            R²ₙ = ε * (norm(v) ^ 2)
            R² = [R²ₙ for k = 1:n]
        else
            return v
        end
    end
end

# ╔═╡ f1cafaab-4605-44f5-a55b-129813277b28
let
    B = [
        63 74 93 93 33
        -14 -20 -46 11 -93
        -1 23 -19 13 12
        84 -32 0 60 57
        61 -52 -63 52 -2
    ]
    @test find_svp_by_enum(B) in ([0, 1, 1, 0, 1], [0, -1, -1, 0, -1])
end

# ╔═╡ b61baecc-19a8-4856-9eac-a14d05c35182
let
    B = [
        -79 43 -1 -58 84 -1 19 -58 17 93
        35 -64 -97 -38 -61 34 16 -17 31 -6
        31 -37 -91 87 93 58 52 99 78 -7
        83 -31 -43 42 -67 -38 32 93 53 -12
        -66 -27 19 94 3 29 -20 -49 40 79
        35 -7 -21 -83 94 67 55 -53 -22 -40
        -32 -42 -65 66 31 -18 94 24 -39 27
        46 21 -36 -69 27 15 -34 51 7 -95
        21 16 34 -2 -60 -75 4 5 70 98
        2 16 -55 -30 98 -16 80 93 -98 20
    ]
    @test find_svp_by_enum(B) in
          ([0, 8, -3, 4, 24, -9, -16, -5, -9, -13], [0, -8, 3, -4, -24, 9, 16, 5, 9, 13])
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
OffsetArrays = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
OffsetArrays = "~1.17.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.1"
manifest_format = "2.0"
project_hash = "ce995e5c35e8d6524e2067c5cd1237b720050f47"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

    [deps.OffsetArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"
"""

# ╔═╡ Cell order:
# ╟─bb3554a0-27c4-4579-a0c6-863a23934709
# ╠═67320cd6-a8a2-11f0-81d5-ff998c22343b
# ╠═851d1600-74b6-4049-8a8d-94638b688e93
# ╠═17033c43-3f58-4cae-9af1-83d87249019d
# ╠═0578cc9c-6620-41b3-b9e3-5304b04eea6b
# ╠═f1cafaab-4605-44f5-a55b-129813277b28
# ╠═b61baecc-19a8-4856-9eac-a14d05c35182
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
