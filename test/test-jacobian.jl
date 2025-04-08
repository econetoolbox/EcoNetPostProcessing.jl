@testset "Jacobian: two producers." begin
    # Default r and K.
    fw = Foodweb([0 0; 0 0])
    m = default_model(fw)
    j = jacobian(m, m.K)
    @test j == [
        -1 0
        0 -1
    ]
    # Non-default r.
    m = default_model(fw, GrowthRate([2, 3]))
    j = jacobian(m, m.K)
    @test j == [
        -2 0
        0 -3
    ]
    # Non-default K.
    m = default_model(fw, CarryingCapacity([2, 3]))
    j = jacobian(m, m.K)
    @test j == [
        -1 0
        0 -1
    ]
    # Producer competition.
    m = default_model(fw, ProducersCompetition(; offdiag = 0.5))
    Beq = [2 / 3, 2 / 3] # Derived analytically.
    j = jacobian(m, Beq)
    @test j == [
        -2/3 -1/3
        -1/3 -2/3
    ]
end

@testset "Jacobian: one producer, one consumer." begin
    # Bioenergetic response.
    fw = Foodweb([0 0; 1 0])
    m = default_model(fw)
    sol = simulate(m, rand(2), 100_000)
    Beq = sol.u[end]
    j = jacobian(m, Beq)
    x, y, h, B0, e = m.x[2], m.y[2], m.h, 0.5, m.e[2, 1]
    F21 = Beq[1]^h / (B0^h + Beq[1]^h)
    dF21dB1 = (h * B0^h * Beq[1]^(h - 1)) / (B0^h + Beq[1]^h)^2
    j21 = x * y * Beq[2] * dF21dB1
    j12 = -F21 * x * y / e
    j11 = (1 - 2Beq[1]) - x * y / e * Beq[2] * dF21dB1
    j22 = 0
    expected = [
        j11 j12
        j21 j22
    ]
    @test isapprox(j, expected, atol = 1e-3)
    # Classic response.
    m = default_model(fw, ClassicResponse())
    sol = simulate(m, rand(2), 100_000)
    Beq = sol.u[end]
    j = jacobian(m, Beq)
    M, a, ht, h = m.M[2], m.attack_rate[2, 1], m.handling_time[2, 1], m.h
    e = m.e[2, 1]
    F21 = a * Beq[1]^h / (M * (1 + a * ht * Beq[1]^h))
    dF21dB1 = (h * a * Beq[1]^(h - 1)) / (M * (1 + a * ht * Beq[1]^h)^2)
    j21 = dF21dB1 * e * Beq[2]
    j22 = 0
    j12 = -F21
    j11 = (1 - 2Beq[1]) - Beq[2] * dF21dB1
    expected = [
        j11 j12
        j21 j22
    ]
    @test isapprox(j, expected, atol = 1e-4)
end

@testset "Jacobian: non-trophic interactions." begin
    # Facilitation.
    fw = Foodweb([0 0; 0 0])
    m = default_model(
        fw,
        ClassicResponse(),
        NontrophicLayers(; :facilitation => (; intensity = 0.1, C = 1)),
    )
    sol = simulate(m, rand(2), 100_000)
    Beq = sol.u[end]
    j = jacobian(m, Beq)
    @test j ≈ [
        -1.1 0
        0 -1.1
    ]
    # Competition.
    fw = Foodweb([0 0; 0 0])
    c0 = 0.1
    m = default_model(
        fw,
        ClassicResponse(),
        NontrophicLayers(; :competition => (; intensity = c0, C = 1)),
    )
    sol = simulate(m, rand(2), 100_000)
    Beq = sol.u[end]
    j = jacobian(m, [0, 2])
    @test j ≈ [
        1-2c0 0
        0 -3
    ]
    # Refuge.
    fw = Foodweb([0 0 0; 0 0 0; 1 0 0])
    r0 = 0
    m = default_model(
        fw,
        ClassicResponse(),
        NontrophicLayers(; :refuge => (; intensity = r0, C = 1)),
    )
    j_norefuge = jacobian(m, [1, 1, 1])
    r0 = 0.1
    m = default_model(
        fw,
        ClassicResponse(),
        NontrophicLayers(; :refuge => (; intensity = r0, C = 1)),
    )
    j_refuge = jacobian(m, [1, 1, 1])
    @test j_refuge[1, 2] > 0
    @test j_norefuge[1, 2] == 0
    # Interspecific interference.
    fw = Foodweb([0 0 0; 1 0 0; 1 0 0])
    i0 = 0.0
    m = default_model(
        fw,
        ClassicResponse(),
        NontrophicLayers(; :interference => (; intensity = i0, C = 1)),
    )
    sol = simulate(m, rand(3), 100_000)
    j_nointerference = jacobian(m, [1, 1, 1])
    i0 = 0.1
    m = default_model(
        fw,
        ClassicResponse(),
        NontrophicLayers(; :interference => (; intensity = i0, C = 1)),
    )
    sol = simulate(m, rand(3), 100_000)
    j_interference = jacobian(m, [1, 1, 1])
    @test j_nointerference[2, 3] == j_nointerference[3, 2] == 0
    @test j_interference[2, 3] == j_interference[3, 2] < 0
end

@testset "Jacobian: nutrients" begin
    fw = Foodweb([0 0; 0 0])
    nutrients = NutrientIntake(;
        half_saturation = [0.3 0.9; 0.9 0.3],
        turnover = 0.9,
        supply = 10,
        concentration = 1,
        r = [1, 2], # Plant intrinsic growth rates.
    )
    m = default_model(fw, nutrients, Mortality([0.6, 1.2]))
    sol = simulate(m, [1, 1], 10_000; N0 = [1, 1])
    j = jacobian(m, sol.u[end])
    @test all(real.(eigvals(j)) .< 1e-3)
end

@testset "Resilience" begin
    fw = Foodweb([0 0; 0 0])
    m = default_model(fw)
    j = jacobian(m, m.K)
    @test resilience(j) ≈ -1
    @test resilience([-0.1 0; 0 -1]) == -0.1
    @test resilience([-0.2 2; 0 -1]) == -0.2
end

@testset "Resilience" begin
    fw = Foodweb([0 0; 0 0])
    m = default_model(fw)
    j = jacobian(m, m.K)
    @test reactivity(j) ≈ -1
    @test reactivity([-0.1 0; 0 -1]) == -0.1
    @test reactivity([-1 -4; 0 -2]) > 0 # Reactive system.
end
