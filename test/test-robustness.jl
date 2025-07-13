@testset "Secondary extinctions" begin
    # Trophic chain - two species.
    fw = Foodweb([2 => 1]) # 2 eats 1.
    m = default_model(fw)
    sol = simulate(m, [1, 1], 1_000)
    B0 = sol.u[end]
    @test secondary_extinctions(m, 1, B0) == [2]
    @test secondary_extinctions(m, 2, B0) == []
    # Trophic chain - three species.
    fw = Foodweb([3 => 2, 2 => 1])
    m = default_model(fw)
    sol = simulate(m, [1, 1, 1], 1_000)
    B0 = sol.u[end]
    @test secondary_extinctions(m, 1, B0) == [2, 3]
    @test secondary_extinctions(m, 2, B0) == [3]
    @test secondary_extinctions(m, 3, B0) == []
end

@testset "Robustness" begin
    n = rand(2:10)
    fw = Foodweb(zeros(Int, n, n))
    m = default_model(fw)
    @test robustness(m) == 0
    fw = Foodweb([2 => 1])
    m = default_model(fw)
    @test 0.4 <= robustness(m) <= 0.6
end
