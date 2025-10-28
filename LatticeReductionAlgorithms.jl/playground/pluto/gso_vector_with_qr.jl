### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 91cf31ee-a3ed-11f0-2f3c-dfe4e6939bb4
using LinearAlgebra

# ╔═╡ 5d08efce-a966-40c8-88d7-760b0197ca0c
md"""
# GSO ベクトルと GSO 係数の計算

GSO は正規化しないグラムシュミットの直交化法による QR 分解を行うことに対応している．
GSO 係数は QR 分解をして得られた上三角行列 $R$ の各行の対角成分で各行を割った行列 $\tilde{R}$ で実現できる．
GSO ベクトルは QR 分解をして得られた $Q$ に対して $R$ の対角成分からなる行列 $D$ を右からかけた行列 $\tilde{Q}$ で実現できる．
"""

# ╔═╡ c38c7e40-2fce-40e2-bd1d-49ba62a5ac34
B = [
	5  2  3
	-3 -7 -10
	-7 -7  0
]

# ╔═╡ 80b29f37-eb9a-4b91-81bc-57bf1bd8bdb9
begin
	F = qr(B)
	R̃ = F.R
	for i in axes(R̃, 1)
		dᵢ = R̃[i, i]
		for j in i:size(R̃, 2)
			R̃[i, j] /= dᵢ
		end
	end
	R̃
end

# ╔═╡ 2fce077d-1094-4a7c-83b8-1d3b7c91501c
Q̃ = F.Q * Diagonal(F.R)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.7"
manifest_format = "2.0"
project_hash = "ac1187e548c6ab173ac57d4e72da1620216bce54"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╟─5d08efce-a966-40c8-88d7-760b0197ca0c
# ╠═91cf31ee-a3ed-11f0-2f3c-dfe4e6939bb4
# ╠═c38c7e40-2fce-40e2-bd1d-49ba62a5ac34
# ╠═80b29f37-eb9a-4b91-81bc-57bf1bd8bdb9
# ╠═2fce077d-1094-4a7c-83b8-1d3b7c91501c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
