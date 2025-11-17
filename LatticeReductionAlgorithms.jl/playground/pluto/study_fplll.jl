### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ a767985f-0022-4c35-8a54-14ac94485336
using DelimitedFiles

# ╔═╡ 8bee8e64-4d46-410c-a505-6e03c829deca
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

# ╔═╡ d93c4eb9-bad7-486b-8b92-7e4fb866b4d8
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
    mat
end

# ╔═╡ c77dd9eb-5984-46f2-97fb-ff35da341e64
begin
    d = 10
    b = 10
end

# ╔═╡ aa37f260-99e2-489a-954d-aa9f34737e43
B = latticegen(d, b)

# ╔═╡ 833e4bef-a637-48c5-ac99-297a969914c5
fplll(d, b)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.1"
manifest_format = "2.0"
project_hash = "c2bd3517fa16afe10381a499d946223d1765af60"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"
"""

# ╔═╡ Cell order:
# ╠═a767985f-0022-4c35-8a54-14ac94485336
# ╠═8bee8e64-4d46-410c-a505-6e03c829deca
# ╠═d93c4eb9-bad7-486b-8b92-7e4fb866b4d8
# ╠═c77dd9eb-5984-46f2-97fb-ff35da341e64
# ╠═aa37f260-99e2-489a-954d-aa9f34737e43
# ╠═833e4bef-a637-48c5-ac99-297a969914c5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
