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
	    # Perform QR on a BigFloat copy with sufficient precision to avoid Inf/NaN
	    # caused by converting very large BigInts to Float64.
	    maxbits = mapreduce(x -> ndigits(abs(x); base=2), max, B; init=0)
	    prec_bits = max(256, 2 * maxbits + 32)
	    # Uniformly scale B by a power-of-two so QR stays in a safe exponent range.
	    # This preserves μ and all ratios used by the algorithms.
	    scale_pow = max(0, maxbits - 256)
	    Ffac = setprecision(prec_bits) do
	        scale = ldexp(BigFloat(1), -scale_pow)
	        qr(scale .* BigFloat.(B))
	    end
	    Rfact = Matrix(Ffac.R)
	    if any(x -> !isfinite(x), Rfact)
	        println("[debug] R from QR contains non-finite entries; sample R[1,2]=", Rfact[1, min(2, size(Rfact,2))])
	        println("[debug] R diag = ", collect(diag(Rfact)))
	    end
	    # Preserve original diagonal for Gram-Schmidt lengths and Q scaling
	    d = diag(Rfact)
	    # Build μ (size-reduction coefficients) safely without in-place mutation
	    μ = copy(Rfact)
	    for i in axes(μ, 1)
	        dᵢ = d[i]
	        if dᵢ == 0 || dᵢ != dᵢ # guard zero/NaN
	            # Zero-norm GS vector: set row to represent no coupling
	            for j = i:size(μ, 2)
	                μ[i, j] = (j == i) ? one(eltype(μ)) : zero(eltype(μ))
	            end
	        else
	            for j = i:size(μ, 2)
	                μ[i, j] = μ[i, j] / dᵢ
	            end
	        end
	    end
	    # Construct GS vectors: Q * diag(original R diagonals)
	    Q̃ = Ffac.Q * Diagonal(d)
	    B⃗ = [dot((@view Q̃[:, j]), (@view Q̃[:, j])) for j in axes(Q̃, 2)]
	    return GSOData(Matrix(B), B⃗, Q̃, μ)
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

	@inline function iszerovec(v)
		r = true
		for i in eachindex(v)
			r *= v[i] ≈ zero(v[i])
		end
		r
	end
	
	function partial_size_reduce!(B::AbstractMatrix{T}, R::AbstractMatrix, i::Int, k::Int) where {T}
	    (1 ≤ i < k) || return
	    μ = R[i, k]
	    m = round(T, μ)
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
function BKZ_reduction!(B::AbstractMatrix, β::Integer, δ::Real)
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
		if norm(g.Q[:, k]) > sqrt(πₖv_norm2) + eps()
			z = 0
			Bsub = Matrix{BigInt}(undef, size(B, 1), h+1)
			for i in 1:k-1
				Bsub[:, i] .= @view B[:, i]
			end
			Bsub[:, k] = v
			for i in k:h
				Bsub[:, 1+i] = @view B[:, i]
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

# ╔═╡ 1f2a9445-9a43-4490-9ddc-eab6a7f86b1e
let
	B = [
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
	volL = abs(det(B))
	@show volL
	n = size(B, 1)
	minkowski_upper_bound = √n * (volL) ^ (1/n)
	@show minkowski_upper_bound
	BKZ_reduction!(B, 10, 0.5)
	minkowski_upper_bound
	@show norm(B[:, 1])
end

# ╔═╡ d5ffe339-8d7f-4fc9-a581-d9491e99b886
let
	B = [
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
	volL = abs(det(B))
	n = size(B, 1)
	minkowski_upper_bound = √n * (volL) ^ (1/n)
	@code_warntype BKZ_reduction!(B, 10, 0.5)
end

# ╔═╡ 1bcdee95-e7d6-4118-a4b0-69873ff1c491
# ╠═╡ disabled = true
#=╠═╡
let
	B = [[2116403082371869720683693394276970360299642162558342355675014410771989791845340105007240644994094893408332181197956887657 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[371644562438531748585630085667746983470408244373580742281269821738816154869236300428137729771977309321210236491639933332 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[970128588937842140661210199069589679615703684562655262446230094731888177573698872897660975186881893815551083050469479861 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[1617885213873858163065250499756376305837976609435400008057297118650398723827143214937143752380433141030895046085134663139 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[861650440856876922136315403300972806928844583476323346957866774048720276102084753290592827764541119558381900492846561687 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[1238401879408921089844474386903737818276239209740348854610933282304058072380604156331847453869318599590680629601714021166 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[302328193892196204243023116106370462383132199983296195152211088133814423239197999286069616594274614510736145263925517783 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[443681408734567838343607093614379754005413471877157728597635239760454843085205167872415159913459109353288267514781828878 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[766924867740885236236436059448427698090379772276332984401846634621581522685296543079817924330812470098383411173196131706 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[1297358997739527278188204320802855592592899794686777403497982861534586993960747704855480346961487202646587358008424555998 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[1923256447555821027760823712014018288715572436486045654757604016247481734877091777864720799782922934671078042956616750575 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[138791684593421316484912136928975449358519306209923545534576080943933941788220318189376736295439954669622762228268690319 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[1752053685015713567891789191640042337741805857685242642164798009587123835789029180360779591122248108173121788944806990018 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[2046303652643607551083773382245919660376167541225343779851326688493989336785866602067814167781776230842039502024961224990 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[161701314003903527294539093851830420185672588693938679547737547227941952578846685981290488906797666486385636306427956922 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[2059714413420462412604210521069862437228594984719691499291824764864200097906908418575923366701593028750494787160369423664 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[236619910841436753359639417128223769093980125208659184950079289394351446700013212860226717401403083201476790108710198995 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[521807735461695034807439230274452751978618453371762238147892673540779228084610103951896314431773738642862361617950116793 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[1898466106725110481310768025369914138166711138228177956551055744552328440573249624614143156799333166991080740657486025033 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[2043257177335956951070219397858777123559357103945545237417070888207682821645422085875324039671029805011652245597869305538 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[904225163671074382768019753539274681795779625382511961264465653291127670514453520050551644296901230356912243654225654401 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[32771872907828211382289261014117057923327418144915569087238809343521061214616361304729214156019094728432041297468304493 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[185789413730166201171458176468755969088909804647627433894406398430754979813498091666025324336151794710682026172503101451 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[832962803741579150035009987801883256172874768041791705764900197393212435979706525525238578679210561346645899576161141793 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[1042621519072433920283083562428469593602182919945835211322354565330887090220897472063018488024892773750305996602609483380 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[171388995082154233233012548090942599143715724202696890942088990082340274801161970138190639906868753248272491463187401695 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[2100493671200104536533474921602153259773348308608862002579228988953071977809598351863259953863461834410595513598037485038 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0]
[849028874593296703934398545250196283716795306401228544442045159484861904205470885973430411009159267090428932112495507852 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]
[1046701824957941065045627030825750373570344803958334982369286898950041088016988929016615138537061510672807965422152953094 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]
[542875223055107687983299240475304624160994620138110511564200316784444891282138558584574906974554732933893987639864801905 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]
[1584382682441666539891259939997296576922760554674178579888435692973226441994929350838648607785926118636355860195320578386 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]
[316503555252451563889967315570982477067169012390726271667771540315451791264064871351402596704077496942199692552980634552 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0]
[1837738370141979294787389742951594377674155679679663044021277108096701623188534975556948353238304199946688862260850711361 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]
[739578827111769709363025786954762546314747907261948920241963639022818938531742051897853101539916787244178755946187647555 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0]
[1619290331242112257658768582175275050536918255267445777163956250682374601071625719555758108294380161274992932056486622464 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0]
[852244520379043286720418573598354118899465075623425806989984806937105488837636350325384144744195783077780523735920132437 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]
[491483555886899765880626820064070785891162464132735054675483057099729053163156386346862511333336434543264870015359686516 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]
[842211543902380452313476978966614900223458281660135345978369910918684855621079286840261544398613444576199814061553688778 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]
[1874411128584029818874533765431554481569149176662206972692659544113727361200950942823993552649265946057735986943133343384 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]
[1807564784681621941238943795385610371403525911952667414563776455384585295644151678837176854600997116364524617338835884611 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]
] |> transpose |> Matrix

	volL = abs(det(B))
	@show volL
	n = size(B, 1)
	minkowski_upper_bound = √n * (volL) ^ (1/n)
	@show minkowski_upper_bound
	BKZ_reduction!(B, 10, 0.5)
	minkowski_upper_bound
	@show norm(B[:, 1])
end
  ╠═╡ =#

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
# ╠═1f2a9445-9a43-4490-9ddc-eab6a7f86b1e
# ╠═d5ffe339-8d7f-4fc9-a581-d9491e99b886
# ╠═1bcdee95-e7d6-4118-a4b0-69873ff1c491
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
