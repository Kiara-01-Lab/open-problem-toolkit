### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 96ab923a-afb1-11f0-86a8-9361eb87280f
begin
	using Test
	using OffsetArrays: OffsetVector
	using LinearAlgebra
end

# ╔═╡ 664f5e2e-ee6e-417d-bbef-2ebd485f0708
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
		return GSOData(Matrix(B), B⃗, Q̃, R̃)
	end
	
	function partial_size_reduce!(g::GSOData{IType, FType}, i::Int, j::Int) where {IType, FType}
		if !(i < j)
			throw(
				ArgumentError("should satisfy i < j, actual i=$(i), j=$(j)")
			)
		end
		μ_ij = g.R[i, j]
		q = fplll_round(μ_ij)
		bi = @view g.B[:, i]
		bj = @view g.B[:, j]
		@. bj -= q * bi
		for l in 1:i
			g.R[l, j] -= q * g.R[l, i]
		end
		g
	end

	@inline function iszerovec(v)
		r = true
		for i in eachindex(v)
			r *= v[i] ≈ zero(v[i])
		end
		r
	end

	function fplll_round(μ)
	    # fplll の最近整数関数に合わせる: 最近接（0.5 はゼロから遠い側に丸め）
	    if μ >= zero(μ)
	        q_i = floor(BigInt, μ + (one(μ) / 2))
	    else
	        q_i = ceil(BigInt, μ - (one(μ) / 2))
	    end
	end

	
	function partial_size_reduce!(B::AbstractMatrix{T}, R::AbstractMatrix, i::Int, k::Int) where {T}
	    (1 ≤ i < k) || return
	    μ = R[i, k]
	    m = fplll_round(μ)
	    m == 0 && return
	    B[:, k] .-= m .* B[:, i]
	    @inbounds begin
	        for j = 1:i-1
	            R[j, k] -= m * R[j, i]
	        end
	        R[i, k] -= m
	    end
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

# ╔═╡ 9acc4eb8-bbce-445a-8775-e34d8dfe9d4b
function ENUM_reduce(
	μ::AbstractMatrix{T}, 
	B⃗::AbstractVector{T}, 
	R²::AbstractVector
)::Tuple{Vector{BigInt}, Bool} where {T <: AbstractFloat}
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
			v[k] = fplll_round(c[k])
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

# ╔═╡ 86981bab-8291-4406-87c2-d0b997a9de90
begin
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
end

# ╔═╡ c6eb0aa3-c592-43dd-96cd-db86c9cc157a
begin
	function _MLLL_reduce!(ℬ::AbstractMatrix{T}, δ::Float64) where {T}
	    h = size(ℬ, 2)
	    z = h
	    g = 1
	
	    # GSO 用のワーク
	    Q = float(T).(ℬ)
	    R = zeros(eltype(Q), size(ℬ))
	    for i in 1:min(size(R, 1), size(R, 2))
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
	        for i = 1:g-1
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
	
	            ν = R[k - 1, k]
	            B = B⃗[k] + ν^2 * B⃗[k - 1]
	
	            if B ≥ δ * B⃗[k - 1]
	                for j = k - 2:-1:1
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
	                    ℬ[:, k - 1], ℬ[:, k] = ℬ[:, k], ℬ[:, k - 1]
	                    for j = 1:k-2
	                        R[j, k], R[j, k - 1] = R[j, k - 1], R[j, k]
	                    end
	
	                    if !(B ≈ 0)
	                        if B⃗[k] ≈ 0
	                            B⃗[k] = B
	                            Q[:, k - 1] .= ν .* Q[:, k - 1]
	                            R[k - 1, k] = inv(ν)
	                            for i = k+1:l
	                                R[k - 1, i] /= ν
	                            end
	                        else
	                            t = B⃗[k - 1] / B
	                            R[k - 1, k] = ν * t
	                            w = Q[:, k - 1]
	                            Q[:, k - 1] .= Q[:, k] .+ ν .* w
	                            B⃗[k - 1] = B
	                            if k ≤ l
	                                Q[:, k] .= -R[k - 1, k] .* Q[:, k] .+ (B⃗[k] / B) .* w
	                                B⃗[k] *= t
	                            end
	                            for i = k+1:l
	                                t = R[k, i]
	                                R[k, i]     = R[k - 1, i] - ν * t
	                                R[k - 1, i] = t + R[k - 1, k] * R[k, i]
	                            end
	                        end
	                    else
	                        B⃗[k], B⃗[k - 1] = B⃗[k - 1], B⃗[k]
	                        Q[:, k], Q[:, k - 1] = Q[:, k - 1], Q[:, k]
	                        for i = k+1:l
	                            R[k, i], R[k - 1, i] = R[k-1, i], R[k, i]
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
			for j in size(B, 2):-1:1
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
end

# ╔═╡ a6c717b5-0ba2-4fa4-9da7-2ed8950e8cef
function find_svp_by_enum(B, k, l)
	@assert k < l
	ε = 0.99
	g = GSOData(B)
	n = size(B, 2)

	R²ₙ = ε * maximum(g.B⃗[k:l])
	R² = [R²ₙ for k in k:l]
	
	μ = g.R[k:l,k:l]
	B⃗ = g.B⃗[k:l]
	coeff, is_succeeded = ENUM_reduce(μ, B⃗, R²)
	return coeff
end

# ╔═╡ fad8d91f-0216-4712-ac84-cb15cf5cb905
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
			@info "coeff got zero vector"
			break
		end
		v = zeros(eltype(B), n)
		for (i, idx_k_to_l) in enumerate(k:l)
			v += coeff[i] * B[:, idx_k_to_l]
		end

		g = GSOData(B)
		
		πₖv_norm2 = 0
		for i = k:n
			πₖv_norm2 += dot(v, g.Q[:, i]) ^ 2 / dot(g.Q[:, i], g.Q[:, i])
		end
		if norm(g.Q[:, k]) > sqrt(πₖv_norm2) + 0.0001
			z = 0
			Bsub = hcat((B[:, i] for i in 1:k-1)..., v, (B[:, i] for i in k:h)...)
			MLLL_reduce!(Bsub, δ)
			B[:, 1:h] .= Bsub[:, 1:h]
		else
			z += 1
			g′ = LLL_reduce((@view B[:, 1:h]), δ)
			B[:, 1:h] = g′.B
		end
	end # while
end

# ╔═╡ 209d5b79-f621-40da-bb56-2cb6f6d479e2
let
	B = [
		63  74  93  93 33
		-14 -20 -46 11 -93
		-1  23  -19 13 12
		84  -32   0 60 57
		61  -52 -63 52 -2
	]
	BKZ_reduce!(B, 3, 0.75)
	B
end 

# ╔═╡ dc1f8ffa-4fca-4d68-8053-a17e7aab3586
let
	B = BigInt[
		    -79   43   -1  -58   84   -1   19  -58   17   93
		     35  -64  -97  -38  -61   34   16  -17   31   -6
		     31  -37  -91   87   93   58   52   99   78   -7
		     83  -31  -43   42  -67  -38   32   93   53  -12
		    -66  -27   19   94    3   29  -20  -49   40   79
		     35   -7  -21  -83   94   67   55  -53  -22  -40
		    -32  -42  -65   66   31  -18   94   24  -39   27
		     46   21  -36  -69   27   15  -34   51    7  -95
		     21   16   34   -2  -60  -75    4    5   70   98
		      2   16  -55  -30   98  -16   80   93  -98   20
	]
	@info minimum([norm(B[:, i]) for i in axes(B, 2)])
	BKZ_reduce!(B, 3, 0.75)
	norm(B[:, 1])
	B[:, 1]
end

# ╔═╡ e22c6910-e1a9-473a-a2f0-fe2aacb9dfb1
let
	B = [
				-2  3  2  8
		 		 7 -2 -8 -9
				 7  6 -9  6
				-5 -1 -7 -4
	]
	δ = 0.9999999
	β = 2
	BKZ_reduce!(B, β, δ)
	B
end

# ╔═╡ 4c47736f-4125-482d-a735-c2492219f652
begin
	function svpchallengedim40seed0()
		str = read(`cat svpchallengedim40seed0.txt`, String)
	    # 1) 先頭 [[ と末尾 ]] を削除
	    s = replace(str, r"^\s*\[\[|\]\n\]\s*$" => "")
	    # 2) ][（または ] [）を行区切り ; に
	    s = replace(s, r"\]\s*\[" => ";")
	    # 3) 1組の角括弧で包む（列はスペース区切りのまま）
	    julia_expr = "[" * s * "]"
	
	    mat = Matrix(transpose(eval(Meta.parse(julia_expr))))
	    mat
	end
	
	function fplll_bkz(β, δ)
		s = read(
			pipeline(
				`cat svpchallengedim40seed0.txt`, 
				`fplll -a bkz -b $(β) -d $(δ)`), 
			String
		)
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

	function fplll_lll(δ)
		s = read(
			pipeline(
				`cat svpchallengedim40seed0.txt`, 
				`fplll -a lll -d $(δ)`), 
			String
		)
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

# ╔═╡ 40583a09-7444-4719-9330-7168d6959887
let
	δ = 0.99
	setprecision(2048) do
		B = svpchallengedim40seed0()
		@time refB = fplll_lll(δ)
		@time LLL_reduce(B, δ)
		@show Float64(norm(refB[:, 1]))
		@show Float64(norm(B[:, 1]))
		T = B \ refB
		T = round.(BigInt, T)
		@assert det(T) in [-1, 1]
	end
end

# ╔═╡ 4c53633b-e60d-437f-88a5-3ada891280d9
let
	d = 40
	b = 100
	β = 20
	δ = 0.99

	setprecision(512) do
		B = svpchallengedim40seed0()
		@time refB = fplll_bkz(β, δ)
		@time BKZ_reduce!(B, β, δ)
		@show Float64(norm(refB[:, 1]))
		@show Float64(norm(B[:, 1]))
		T = B \ refB
		T = round.(BigInt, T)
		@assert det(T) in [-1, 1]
	end
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
# ╠═96ab923a-afb1-11f0-86a8-9361eb87280f
# ╠═664f5e2e-ee6e-417d-bbef-2ebd485f0708
# ╠═9acc4eb8-bbce-445a-8775-e34d8dfe9d4b
# ╠═86981bab-8291-4406-87c2-d0b997a9de90
# ╠═c6eb0aa3-c592-43dd-96cd-db86c9cc157a
# ╠═a6c717b5-0ba2-4fa4-9da7-2ed8950e8cef
# ╠═fad8d91f-0216-4712-ac84-cb15cf5cb905
# ╠═209d5b79-f621-40da-bb56-2cb6f6d479e2
# ╠═dc1f8ffa-4fca-4d68-8053-a17e7aab3586
# ╠═e22c6910-e1a9-473a-a2f0-fe2aacb9dfb1
# ╠═4c47736f-4125-482d-a735-c2492219f652
# ╠═40583a09-7444-4719-9330-7168d6959887
# ╠═4c53633b-e60d-437f-88a5-3ada891280d9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
