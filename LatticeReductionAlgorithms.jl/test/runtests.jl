using Test

using LinearAlgebra
using LatticeReductionAlgorithms: GSOData, LLL_reduce, BKZ_reduce!, MLLL_reduce!

@testset "GSOData" begin
    B = [
        5 2 3
        -3 -7 -10
        -7 -7 0
    ]
    g = GSOData(B)

    (; Q, R) = g
    @testset "直交性" begin
        for i in axes(R, 2)
            for j = (i+1):size(R, 2)
                @test abs(dot(Q[:, i], Q[:, j])) < 1e-13
            end
        end
    end

    @testset "norm" begin
        for i in axes(R, 2)
            @test norm(Q[:, i]) ≤ norm(B[:, i])
        end
    end

    @testset "volume" begin
        volL = abs(det(B))
        @test volL ≈ prod(norm(Q[:, i]) for i in axes(R, 2))
    end
end

@testset "LLL_reuce" begin
    @testset "Example 2.3.9 for δ = 0.75" begin
        B = [
            9 8 3
            2 6 2
            7 1 6
        ]
        δ = 0.75
        g = LLL_reduce(B, δ)
        @test g.B == [
            -1 2 3
            4 6 -2
            -6 0 -5
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
            6 3 2
            0 -2 6
            1 -5 0
        ]
    end

    @testset "Example 2.3.10 for δ = 0.9999999" begin
        B = [
            -2 3 2 8
            7 -2 -8 -9
            7 6 -9 6
            -5 -1 -7 -4
        ]
        δ = 0.9999999
        g = LLL_reduce(B, δ)
        @test g.B == [
            2 2 -2 3
            3 0 2 -2
            1 -2 3 6
            1 -4 -3 -1
        ]

        α = 4 / (4δ - 1)
        n = size(B, 2)
        volL = abs(det(B))
        @test norm(g.B[:, 1]) ≤ α ^ ((n - 1)/4) * volL ^ (1/n)
        @test prod(norm(g.B[:, i]) for i = 1:n) ≤ α ^ (n * (n - 1)/4) * volL
    end
end

@testset "MLL_reduce!" begin
	δ = 0.75
	ℬ = [
		388 -672 -689 -179 508
		417 -73  379  96   -705
		417 -121 724  -24  173
		-86 944  653  978  -343
	]
	MLLL_reduce!(ℬ, δ)
	expected = [
		 -1   1  -1   3  0
		  0  -1  -4  -1  0
		 -1  -1   2  -3  0
		  0  -3   0   2  0
	]
	@test ℬ == expected
	@test det(ℬ[:, 1:minimum(size(ℬ))]) ≠ 0

	δ = 0.75
	ℬ = [
		-696  -760 552 -160 307  117
		-186  -106 6   -439 -526 -94
		661   -775 9   -544 862  472
		-727   659 726  365 396  138
	]
	MLLL_reduce!(ℬ, δ)
	@test ℬ == [
		 1   0   0   0  0  0
		 0   1  -1  -1  0  0
		 0  -1   0  -1  0  0
		 0   0   1  -1  0  0
	]
	@test det(ℬ[:, 1:minimum(size(ℬ))]) ≠ 0
end

@testset "BKZ_reduce!" begin
    B = BigInt[
        63 74 93 93 33
        -14 -20 -46 11 -93
        -1 23 -19 13 12
        84 -32 0 60 57
        61 -52 -63 52 -2
    ]
    BKZ_reduce!(B, 3, 0.75)
    @test B[:, 1] == [0, 1, 1, 0, 1]

    B = BigInt[
        -79 43 -1 -58 84 -1 19 -58 17 93
        35 -64 -97 -38 -61 34 16 -17 31 -6
        31 -37 -91 87 93 58 52 99 78 -7
        83 -31 -43 42 -67 -38 32 93 53 -12
        -66 -27 19 94 3 29 -20 -49 40 79
        35 -7 -21 -83 94 67 55 -53 -22 -40
        -32 -42 -65 66 31 -18 94 24 -39 27
        46 21 -36 -69 27 15 -34 51 7 -95
        21 16 34 -2 -60 -75 4 5 70 98
        2 16 -55 -30 98 -16 80 93 -98 20
    ]
    BKZ_reduce!(B, 3, 0.75)
    @test B[:, 1] == [
        -2
        14
        -8
        -5
        -11
        -9
        -5
        5
        -18
        21
    ]
end
