### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 11cfb05a-aefa-11f0-bb95-6fb5f8337028
begin
    using Test
    using LinearAlgebra
end

# ╔═╡ 9e1cb7bf-f309-406c-8d6f-7f8004bfb729
begin
    @inline function iszerovec(v)
        r = true
        for i in eachindex(v)
            r *= v[i] ≈ zero(v[i])
        end
        r
    end

    function partial_size_reduce!(
        B::AbstractMatrix{T},
        R::AbstractMatrix,
        i::Int,
        k::Int,
    ) where {T}
        (1 ≤ i < k) || return
        μ = R[i, k]
        m = round(T, μ, RoundNearestTiesAway)
        m == 0 && return
        B[:, k] .-= m .* B[:, i]
        @inbounds begin
            for j = 1:(i-1)
                R[j, k] -= m * R[j, i]
            end
            R[i, k] -= m
        end
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
end

# ╔═╡ e400d08f-4173-4534-a256-813842ff7dd8
let
    δ = 0.75
    ℬ = [
        388 -672 -689 -179 508
        417 -73 379 96 -705
        417 -121 724 -24 173
        -86 944 653 978 -343
    ]
    MLLL_reduce!(ℬ, δ)
    expected = [
        -1 1 -1 3 0
        0 -1 -4 -1 0
        -1 -1 2 -3 0
        0 -3 0 2 0
    ]
    @test ℬ == expected
    @test det(ℬ[:, 1:minimum(size(ℬ))]) ≠ 0
end

# ╔═╡ 86e5aa33-3eee-4da5-b10d-969b500c8a3c
let
    δ = 0.75
    ℬ = [
        -696 -760 552 -160 307 117
        -186 -106 6 -439 -526 -94
        661 -775 9 -544 862 472
        -727 659 726 365 396 138
    ]
    MLLL_reduce!(ℬ, δ)
    @test ℬ == [
        1 0 0 0 0 0
        0 1 -1 -1 0 0
        0 -1 0 -1 0 0
        0 0 1 -1 0 0
    ]
    @test det(ℬ[:, 1:minimum(size(ℬ))]) ≠ 0
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
# ╠═11cfb05a-aefa-11f0-bb95-6fb5f8337028
# ╠═9e1cb7bf-f309-406c-8d6f-7f8004bfb729
# ╠═e400d08f-4173-4534-a256-813842ff7dd8
# ╠═86e5aa33-3eee-4da5-b10d-969b500c8a3c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
