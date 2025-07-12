@testset "Secondary extinctions" begin
    fw = Foodweb([2 => 1]) # 2 eats 1.
    m = default_model(fw)
    sol = simulate(m, [1, 1], 1_000)
    B0 = sol.u[end]
    @test secondary_extinctions(m, 1, B0) == [2]
    @test secondary_extinctions(m, 2, B0) == []
end
