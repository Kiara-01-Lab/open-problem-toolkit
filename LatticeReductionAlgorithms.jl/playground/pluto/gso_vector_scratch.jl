### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 91cf31ee-a3ed-11f0-2f3c-dfe4e6939bb4
begin
	using Test
	using LinearAlgebra
end

# ╔═╡ 5d08efce-a966-40c8-88d7-760b0197ca0c
md"""
# GSO ベクトルと GSO 係数の計算

Gram-Schmidt の正規化しない直交化法のアルゴリズムを実装する．
"""

# ╔═╡ c38c7e40-2fce-40e2-bd1d-49ba62a5ac34
B = [
	5  2  3
	-3 -7 -10
	-7 -7  0
]

# ╔═╡ 80b29f37-eb9a-4b91-81bc-57bf1bd8bdb9
begin
	Q = Float64.(B)
	R = zeros(size(B))
	for i in axes(R, 1)
		R[i, i] = 1
	end
	for j in axes(R, 2)
		Q[:, j] .= @view B[:, j]
		for i in 1:(j-1)
			b_i_ast = @view Q[:, i]
			b_j = @view B[:, j]
			μ_ij = dot(b_j, b_i_ast)/dot(b_i_ast, b_i_ast)
			R[i, j] = μ_ij
			Q[:, j] .-= μ_ij .* b_i_ast
		end
	end
end

# ╔═╡ 1e9001bb-3e19-459c-8504-a15f675769cb
display(Q)

# ╔═╡ b0e80cee-287b-471a-9378-85bcfe49769f
display(R)

# ╔═╡ d8f3639b-fa04-46d1-8aa6-30bc6d9bd4f0
md"""
GSO 行列の性質の確認をする
"""

# ╔═╡ 34a536bf-58c6-46f3-8ddb-82d97ef4a4fd
@testset "直交性" begin
	for i = axes(R, 2)
		for j in (i+1):size(R, 2)
			@test abs(dot(Q[:, i], Q[:, j])) < 1e-13
		end
	end
end

# ╔═╡ ab36b712-f030-40dd-b8b2-26020ba12054
@testset "norm" begin
	for i = axes(R, 2)
		@test norm(Q[:, i]) ≤ norm(B[:, i])
	end
end

# ╔═╡ 269661a9-95ce-4151-bd69-d5406986b56d
@testset "volume" begin
	volL = abs(det(B))
	@test volL ≈ prod(norm(Q[:, i]) for i in axes(R, 2))
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

julia_version = "1.12.1"
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
version = "5.15.0+0"
"""

# ╔═╡ Cell order:
# ╟─5d08efce-a966-40c8-88d7-760b0197ca0c
# ╠═91cf31ee-a3ed-11f0-2f3c-dfe4e6939bb4
# ╠═c38c7e40-2fce-40e2-bd1d-49ba62a5ac34
# ╠═80b29f37-eb9a-4b91-81bc-57bf1bd8bdb9
# ╠═1e9001bb-3e19-459c-8504-a15f675769cb
# ╠═b0e80cee-287b-471a-9378-85bcfe49769f
# ╟─d8f3639b-fa04-46d1-8aa6-30bc6d9bd4f0
# ╠═34a536bf-58c6-46f3-8ddb-82d97ef4a4fd
# ╠═ab36b712-f030-40dd-b8b2-26020ba12054
# ╠═269661a9-95ce-4151-bd69-d5406986b56d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
