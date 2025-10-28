### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 97313604-b3ca-11f0-b308-b9c15180c612
begin
	using Test
	using LinearAlgebra
	using MultiFloats
end

# ╔═╡ 5b290410-03a6-43d5-8340-8158216057e4
begin
	struct GSOData{I<:Real, F<:Real}
		B::Matrix{I}
		B⃗::Vector{F}
		Q::Matrix{F}
		R::Matrix{F}
	end

	function GSOData(B::AbstractMatrix)
		F = float(eltype(B))
		Q = F.(B)
		R = zeros(F, size(B))
		for i in axes(R, 1)
			R[i, i] = 1
		end
		for j in axes(R, 2)
			Q[:, j] .= @view B[:, j]
			for i in 1:(j-1)
				b_i_ast = @view Q[:, i]
				b_j = @view B[:, j]
				μ_ij = dot(b_j, b_i_ast) / dot(b_i_ast, b_i_ast)
				R[i, j] = μ_ij
				Q[:, j] .-= μ_ij .* b_i_ast
			end
		end
		B⃗ = [dot((@view Q[:, j]), (@view Q[:, j])) for j in axes(Q, 2)]
		return GSOData(B, B⃗, Q, R)
	end
	
	function partial_size_reduce!(g::GSOData{IType, FType}, i::Int, j::Int) where {IType, FType}
		if !(i < j)
			throw(
				ArgumentError("should satisfy i < j, actual i=$(i), j=$(j)")
			)
		end
		μ_ij = g.R[i, j]
		q_ = round(BigInt, BigFloat(μ_ij))
		q = FType(q_)
		bi = @view g.B[:, i]
		bj = @view g.B[:, j]
		@. bj -= q * bi
		for l in 1:i
			g.R[l, j] -= q * g.R[l, i]
		end
		g
	end
	
	function size_reduce!(g::GSOData)
		R = g.R
		for j in 2:size(R, 2)
			for i in (j-1):-1:1
				partial_size_reduce!(g, i, j)
			end
		end
		g
	end

	@testset "GSO vector is invariant" begin
		B = [
			5  2  3
			-3 -7 -10
			-7 -7  0
		]
		Q = GSOData(B).Q
		B_reduced = [
			5  -3 1
			-3 -4 -3
			-7 0 7
		]
		Q_reduced = GSOData(B_reduced).Q
		@test isapprox(Q, Q_reduced, atol=1e-14)
	end
end

# ╔═╡ 2fbeef28-29ad-4ff1-9473-41a8fd1fd8d2
function gsoupdate!(g::GSOData, k::Integer)
	k ≥ 2 || throw(ArgumentError("k should satisfy k ≥ 2, actual k=$(k)"))
	for i in axes(g.B, 1)
		# swap
		g.B[i, k-1], g.B[i, k] = g.B[i, k], g.B[i, k-1]
	end
	μ = g.R
	ν = μ[k-1, k]
	B = g.B⃗[k] + ν ^ 2 * g.B⃗[k-1]
	μ[k-1, k] = ν * g.B⃗[k - 1] / B
	g.B⃗[k] = g.B⃗[k] * g.B⃗[k-1] / B
	g.B⃗[k-1] = B
	for j = 1:(k-2)
		# swap
		μ[j, k-1], μ[j, k] =μ[j, k], μ[j, k-1]
	end
	n = size(μ, 2)
	for i = (k+1):n
		t = μ[k, i]
		μ[k, i] = μ[k-1, i] - ν * t
		μ[k-1, i] = t + μ[k-1, k] * μ[k, i]
	end
	g
end

# ╔═╡ 746e0d8d-ceba-4b9d-9b8a-58ef90f301ba
⪆(a, b) = a ≥ b || a ≈ b

# ╔═╡ 15f48411-c9dc-4b40-a04b-abb134f41cc0
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
		if g.B⃗[k] ⪆ (δ - g.R[k-1, k] ^ 2) * g.B⃗[k-1]
			# Lovász 条件を満たす
			k += 1
		else
			# swap basis
			gsoupdate!(g, k)
			k = max(k-1, 2)
		end
	end
	g
end

# ╔═╡ 6295c9c1-9203-4263-8a22-925f3740fbe8
begin
	function latticegen(d, b)
		str = read(`latticegen u $(d) $(b)`, String)
	    # 1) 先頭 [[ と末尾 ]] を削除
	    s = replace(str, r"^\s*\[\[|\]\]\s*$" => "")
	    # 2) ][（または ] [）を行区切り ; に
	    s = replace(s, r"\]\s*\[" => ";")
	    # 3) 1組の角括弧で包む（列はスペース区切りのまま）
	    julia_expr = "[" * s * "]"
	
	    mat = Matrix(transpose(eval(Meta.parse(julia_expr))))
	    mat
	end
	
	function fplll(d, b)
		s = read(pipeline(`latticegen u $(d) $(b)`, `fplll -a lll -d 0.99`), String)
	    # 1) 文字列中の "\n" を実際の改行に
	    s = replace(s, "\\n" => "\n")
	
	    # 2) ][（ ] [ 含む）を行区切り ; に
	    s = replace(s, r"\]\s*\[" => " ; ")
	
	    # 3) 残る全ての [ と ] を削除（余計な ] があっても消える）
	    s = replace(s, ['[', ']'] => ' ')
	
	    # 4) 余分な空白を正規化
	    s = replace(s, r"[ \t]+" => " ")
	    s = strip(s)
	
	    # 5) 1 組の角括弧で包んで Julia の行列リテラルに
	    expr = "[" * s * "]"
	
	    # デバッグ用に確認したい場合:
	    # println(expr)
	
	    mat = Matrix(transpose(eval(Meta.parse(expr))))
	end
end

# ╔═╡ 98765d50-9c42-45f0-8896-d8fef005d03e
begin
	d = 17
	b = 4
	B = latticegen(d, b)
	g = LLL_reduce(B, 0.99)
	display(g.B - fplll(d, b))
	display(g.B)
	display(fplll(d, b))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
MultiFloats = "bdf0d083-296b-4888-a5b6-7498122e68a5"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
MultiFloats = "~2.3.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.1"
manifest_format = "2.0"
project_hash = "c60856abf658f1600d6fed8b76b3a7fdb2341a9d"

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

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

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

[[deps.MultiFloats]]
deps = ["LinearAlgebra", "Printf", "Random", "SIMD"]
git-tree-sha1 = "39ffa6286f40544ecea725d8031c615e79d88d45"
uuid = "bdf0d083-296b-4888-a5b6-7498122e68a5"
version = "2.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"
"""

# ╔═╡ Cell order:
# ╠═97313604-b3ca-11f0-b308-b9c15180c612
# ╠═5b290410-03a6-43d5-8340-8158216057e4
# ╠═2fbeef28-29ad-4ff1-9473-41a8fd1fd8d2
# ╠═746e0d8d-ceba-4b9d-9b8a-58ef90f301ba
# ╠═15f48411-c9dc-4b40-a04b-abb134f41cc0
# ╠═6295c9c1-9203-4263-8a22-925f3740fbe8
# ╠═98765d50-9c42-45f0-8896-d8fef005d03e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
