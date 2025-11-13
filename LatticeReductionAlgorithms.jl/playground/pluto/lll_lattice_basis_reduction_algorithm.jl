### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ eabf46c8-b098-4bb7-8126-454ee3b87e87
begin
	using Test
	using LinearAlgebra
end

# ╔═╡ 52f0adc8-a4b4-11f0-9f15-f38c58262cb5
md"""
LLL 基底簡約アルゴリズムの実装
"""

# ╔═╡ ebb6f39f-3594-4ae6-b86b-238087b3101a
begin
	struct GSOData{I<:Integer, F<:Real}
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
			for j in i:size(R̃, 2)
				R̃[i, j] /= dᵢ
			end
		end
		Q̃ = F.Q * Diagonal(F.R)
		B⃗ = [dot((@view Q̃[:, j]), (@view Q̃[:, j])) for j in axes(Q̃, 2)]
		return GSOData(B, B⃗, Q̃, R̃)
	end
	
	function partial_size_reduce!(g::GSOData{IType, FType}, i::Int, j::Int) where {IType, FType}
		if !(i < j)
			throw(
				ArgumentError("should satisfy i < j, actual i=$(i), j=$(j)")
			)
		end
		μ_ij = g.R[i, j]
		q = round(IType, μ_ij)
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

# ╔═╡ 21495c6b-b02f-42f0-8c27-c763fefc1711
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

# ╔═╡ cf7f8041-26a5-4e9d-aabe-0b69782c446b
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
	g
end

# ╔═╡ af41f90e-f501-4db7-995f-6e94da1f7725
@testset "LLL_reduce" begin
	@testset "Example 2.3.9 for δ = 0.75" begin
		B = [
			9 8 3
			2 6 2
			7 1 6
		]
		δ = 0.75
		g = LLL_reduce(B, δ)
		@test g.B == [
			 -1  2   3
			  4  6  -2
			 -6  0  -5
		]
	end
	
	@testset "Example 2.3.9 for δ = 0.99" begin
		B = [
			9 8 3
			2 6 2
			7 1 6
		]
		δ = 0.99
		g = LLL_reduce(B, δ)
		@test g.B == [
			6 3  2
			0 -2 6
			1 -5 0
		]
	end
	
	@testset "Example 2.3.10 for δ = 0.9999999" begin
		B = [
			-2  3  2  8
	 		 7 -2 -8 -9
			 7  6 -9  6
			-5 -1 -7 -4
		]
		δ = 0.9999999
		g = LLL_reduce(B, δ)
		@test g.B == [
			 2   2  -2   3
			 3   0   2  -2
			 1  -2   3   6
			 1  -4  -3  -1
		]

		α = 4 / (4δ - 1)
		n = size(B, 2)
		volL = abs(det(B))
		@test norm(g.B[:, 1]) ≤ α ^ ((n - 1)/4) * volL ^ (1/n)
		@test prod(norm(g.B[:, i]) for i in 1:n) ≤ α ^ (n * (n - 1)/4) * volL
	end
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
# ╟─52f0adc8-a4b4-11f0-9f15-f38c58262cb5
# ╠═eabf46c8-b098-4bb7-8126-454ee3b87e87
# ╠═ebb6f39f-3594-4ae6-b86b-238087b3101a
# ╠═21495c6b-b02f-42f0-8c27-c763fefc1711
# ╠═cf7f8041-26a5-4e9d-aabe-0b69782c446b
# ╠═af41f90e-f501-4db7-995f-6e94da1f7725
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
