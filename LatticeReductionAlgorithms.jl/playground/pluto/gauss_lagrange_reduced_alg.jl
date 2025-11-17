### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 92d495ae-ee71-4e41-bc09-90da0d29f86b
using LinearAlgebra

# ╔═╡ 9a7c20fc-9a1c-423b-be7a-b2e3986f880f
md"""
# Gauss Lagrange 基底簡約アルゴリズム

青野，安田. 格子暗号解読のための数学的基礎 を参考に Gauss-Lagrange アルゴリズムを実装する．
このアルゴリズムは二次元格子にのみ適用できるアルゴリズムであり．LLL簡約はこのGauss-Lagrange簡約アルゴリズムの一般化である．
"""

# ╔═╡ b3fd66e2-2b64-4e39-b74e-538148e2b84b
struct GaussLagrange end

# ╔═╡ 355d71c0-5c38-454e-bb42-093142125828
function reducebasis!(::Type{GaussLagrange}, B)
    if norm(@view B[:, 1]) > norm(@view B[:, 2])
        # swap cols
        B .= @view B[:, [2, 1]]
    end
    v_tmp = Vector{eltype(B)}(undef, length(@view B[:, 1]))
    b1 = @view B[:, 1]
    b2 = @view B[:, 2]
    while true
        if norm(b1) > norm(b2)
            break
        end
        q = -round(Int, dot(b1, b2) / norm(b1)^2)
        @. v_tmp = b2 + q * b1
        @. b2 = b1
        @. b1 = v_tmp
    end
    # swap cols
    B .= @view B[:, [2, 1]]
    return B
end

# ╔═╡ e7ccc1d8-9a8d-11f0-364f-8bc51dbac3b3
let
    # Example 2.1.7
    b1 = [
        -7
        -4
        -10
    ]

    b2 = [
        9
        5
        12
    ]
    B = hcat(b1, b2)
    reducebasis!(GaussLagrange, B)

    @assert B == [
        1 2
        0 1
        -2 2
    ]

    # ガウス-ラグランジュ簡約基底の性質の確認
    b1 = @view B[:, 1]
    b2 = @view B[:, 2]
    @assert norm(b1) ≤ norm(b2)
    @assert norm(b2) ≤ norm(b1 - b2)
    @assert norm(b2) ≤ norm(b1 + b2)
end

# ╔═╡ d86c591c-6cbc-49fe-9a25-8dcb54aa28dc
let
    # Example 2.1.7
    b1 = [
        230
        -651
        609
        -366
    ]

    b2 = [
        301
        -852
        797
        -479
    ]
    B = hcat(b1, b2)
    reducebasis!(GaussLagrange, B)
    @assert B == [
        -1 2
        -3 -3
        -2 5
        -1 -2
    ]

    # ガウス-ラグランジュ簡約基底の性質の確認
    b1 = @view B[:, 1]
    b2 = @view B[:, 2]
    @assert norm(b1) ≤ norm(b2)
    @assert norm(b2) ≤ norm(b1 - b2)
    @assert norm(b2) ≤ norm(b1 + b2)
end

# ╔═╡ 5c3a953c-1453-489f-b772-5eb84be98c9b
let
    # See 2dlattice_viewer.jl
    b1 = [
        18
        1
    ]

    b2 = [
        8
        1
    ]
    B = hcat(b1, b2)
    reducebasis!(GaussLagrange, B)
    @assert B == [
        2 2
        -1 4
    ]
    # ガウス-ラグランジュ簡約基底の性質の確認
    b1 = @view B[:, 1]
    b2 = @view B[:, 2]
    @assert norm(b1) ≤ norm(b2)
    @assert norm(b2) ≤ norm(b1 - b2)
    @assert norm(b2) ≤ norm(b1 + b2)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.1"
manifest_format = "2.0"
project_hash = "f352ceee806168c8ae38887a01d7bae6ca62470b"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"
"""

# ╔═╡ Cell order:
# ╟─9a7c20fc-9a1c-423b-be7a-b2e3986f880f
# ╠═92d495ae-ee71-4e41-bc09-90da0d29f86b
# ╠═b3fd66e2-2b64-4e39-b74e-538148e2b84b
# ╠═355d71c0-5c38-454e-bb42-093142125828
# ╠═e7ccc1d8-9a8d-11f0-364f-8bc51dbac3b3
# ╠═d86c591c-6cbc-49fe-9a25-8dcb54aa28dc
# ╠═5c3a953c-1453-489f-b772-5eb84be98c9b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
