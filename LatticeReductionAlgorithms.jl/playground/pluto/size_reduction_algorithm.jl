### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 6206c1ac-a4a6-11f0-25cf-17f96340494a
begin
    using Test
    using LinearAlgebra
end

# ╔═╡ a0b75700-73b0-4907-b827-5a884936ea38
md"""
サイズ簡約アルゴリズムの実装
"""

# ╔═╡ 09c9a0f6-61bb-4858-96be-3fd27c0a25ce
struct GSOData
    B::Matrix{Float64}
    Q::Matrix{Float64}
    R::Matrix{Float64}
end

# ╔═╡ fea0c5e0-94d8-4657-877c-3ac1d037e8ee
function gso(B::AbstractMatrix)
    F = qr(B)
    R̃ = F.R
    for i in axes(R̃, 1)
        dᵢ = R̃[i, i]
        for j = i:size(R̃, 2)
            R̃[i, j] /= dᵢ
        end
    end
    Q̃ = F.Q * Diagonal(F.R)
    return GSOData(B, Q̃, R̃)
end

# ╔═╡ 1bc2d62b-28f7-4642-9aa1-4ea39f02cd97
function partial_size_reduce!(d::GSOData, i::Int, j::Int)
    if !(i < j)
        throw(ArgumentError("should satisfy i < j, actual i=$(i), j=$(j)"))
    end
    μ_ij = d.R[i, j]
    q = round(Int, μ_ij)
    bi = @view d.B[:, i]
    bj = @view d.B[:, j]
    @. bj -= q * bi
    for l = 1:i
        d.R[l, j] -= q * d.R[l, i]
    end
    d
end

# ╔═╡ f179a6a5-4e77-4b50-99c6-411166870b37
function size_reduce!(d::GSOData)
    R = d.R
    for j = 2:size(R, 2)
        for i = (j-1):-1:1
            partial_size_reduce!(d, i, j)
        end
    end
    d
end

# ╔═╡ 5a53f73c-c9f8-480e-a646-6e13e5fa162e
@testset "partial_size_reduce!" begin
    B = [
        5 2 3
        -3 -7 -10
        -7 -7 0
    ]
    g = gso(B)
    g = partial_size_reduce!(g, 1, 2)
    @test g.B == [
        5 -3 3
        -3 -4 -10
        -7 -0 0
    ]
    @test g.R ≈ [
        1.0 -0.03614457831325302 0.5421686746987951
        0.0 1.0 1.3107454017424978
        0.0 0.0 1.0
    ]
    g = partial_size_reduce!(g, 2, 3)
    @test g.B == [
        5 -3 6
        -3 -4 -6
        -7 0 0
    ]
    g = partial_size_reduce!(g, 1, 3)
    @test g.B == [
        5 -3 1
        -3 -4 -3
        -7 0 7
    ]
    # size reduced
    for i = 1:size(g.R, 1)
        for j = (i+1):size(g.R, 2)
            @test abs(g.R[i, j]) ≤ 0.5
        end
    end
end

# ╔═╡ 34a2ba4d-eb59-4bbf-8f89-200035028a44
@testset "partial_size_reduce!" begin
    B = [
        5 2 3
        -3 -7 -10
        -7 -7 0
    ]
    g = gso(B)
    g = size_reduce!(g)
    # size reduced
    @test g.B == [
        5 -3 1
        -3 -4 -3
        -7 0 7
    ]
    for i = 1:size(g.R, 1)
        for j = (i+1):size(g.R, 2)
            @test abs(g.R[i, j]) ≤ 0.5
        end
    end
end

# ╔═╡ 0aab64dc-671d-47e9-8a30-9a05653a41b2
@testset "GSO vector is invariant" begin
    B = [
        5 2 3
        -3 -7 -10
        -7 -7 0
    ]
    Q = gso(B).Q
    B_reduced = [
        5 -3 1
        -3 -4 -3
        -7 0 7
    ]
    Q_reduced = gso(B_reduced).Q
    @test isapprox(Q, Q_reduced, atol = 1e-14)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.0"
manifest_format = "2.0"
project_hash = "daa57a5ae35ddc41f805b1b81e5a5b0171a8b179"

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
version = "5.13.1+1"
"""

# ╔═╡ Cell order:
# ╟─a0b75700-73b0-4907-b827-5a884936ea38
# ╠═6206c1ac-a4a6-11f0-25cf-17f96340494a
# ╠═09c9a0f6-61bb-4858-96be-3fd27c0a25ce
# ╠═fea0c5e0-94d8-4657-877c-3ac1d037e8ee
# ╠═1bc2d62b-28f7-4642-9aa1-4ea39f02cd97
# ╠═f179a6a5-4e77-4b50-99c6-411166870b37
# ╠═5a53f73c-c9f8-480e-a646-6e13e5fa162e
# ╠═34a2ba4d-eb59-4bbf-8f89-200035028a44
# ╠═0aab64dc-671d-47e9-8a30-9a05653a41b2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
