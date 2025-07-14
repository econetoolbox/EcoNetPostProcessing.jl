@testset "Interaction matrix - Sensitivity." begin
    # Two producers, no competition.
    fw = Foodweb([0 0; 0 0])
    m = default_model(fw)
    A = get_interaction(m, rand(2)) # rand because A is not density dependent.
    @test A ≈ [-1 0; 0 -1]
    Beq = simulate(m, [1, 1], 1_000; show_degenerated).u[end]
    S = sensitivity(m, Beq)
    @test S ≈ [1 0; 0 1]
    # Two producers with competition.
    A_expected = [-1 -0.2; -0.3 -1]
    comp = ProducersCompetition(-A_expected)
    m = default_model(fw, comp)
    Beq = simulate(m, [1, 1], 1_000; show_degenerated).u[end]
    A = get_interaction(m, rand(2)) # rand because A is not density dependent.
    @test A ≈ A_expected
    S = sensitivity(m, Beq)
    for sp_press in 1:2
        K_press = CarryingCapacity(m.K - 0.1 * ([1, 2] .== sp_press))
        m_press = default_model(fw, comp, K_press)
        B_press = simulate(m_press, Beq, 1_000; show_degenerated).u[end]
        B_expected = Beq .- 0.1 * S[:, sp_press]
        @test B_press ≈ B_expected
    end
    m2 = copy(m)
    # One producer, one consumer - bioenergetic response.
    fw = Foodweb([0 1; 0 0])
    m = default_model(fw)
    sol = simulate(m, [1, 1], 1_000; show_degenerated)
    Beq = sol.u[end]
    j = jacobian(m, Beq)
    A = get_interaction(m, Beq)
    x, y, e, h = m.x[1], m.y[1], m.e[1, 2], m.h
    B0 = m.half_saturation_density[1]
    A11 = 0
    F21 = Beq[2]^h / (B0^h + Beq[2]^h)
    A21 = -x * y * F21 / e / Beq[2]
    A22 = -1 + x * y / e * Beq[1] * (F21 / Beq[2]^2 - B0^h * h / (B0^h + Beq[2]^h)^2)
    A12 = x * y * h * B0^h * Beq[2]^(h - 1) / (B0^h + Beq[2]^h)^2
    @test A ≈ [
        A11 A12
        A21 A22
    ]
    S = sensitivity(m, Beq)
    d = 0.01 * [0.1, 0.5]
    mortality = Mortality(d)
    m_press = default_model(fw, mortality)
    B_press = simulate(m_press, Beq, 1_000).u[end]
    B_expected = Beq .- S * d
    @test B_press ≈ B_expected atol = 1e-3
    # One producer, one consumer - classic response.
    m = default_model(fw, ClassicResponse())
    sol = simulate(m, [1, 1], 1_000; show_degenerated)
    Beq = sol.u[end]
    j = jacobian(m, Beq)
    A = get_interaction(m, Beq)
    @test Diagonal(Beq) * A ≈ j
    S = sensitivity(m, Beq)
    d = 0.01 * [0.1, 0.5]
    mortality = Mortality(d)
    m_press = default_model(fw, mortality, ClassicResponse())
    B_press = simulate(m_press, Beq, 1_000).u[end]
    B_expected = Beq .- S * d
    @test B_press ≈ B_expected atol = 1e-3
end

@testset "Keystoneness." begin
    # Two independent species.
    fw = Foodweb([0 0; 0 0])
    m = default_model(fw)
    @test keystoneness(m, m.K) == [0, 0] # Null for independent species.
    # Two species with competition.
    m = default_model(fw, ProducersCompetition([1 0.1; 0.3 1]))
    Beq = simulate(m, [1, 1], 1_000; show_degenerated).u[end]
    S = sensitivity(m, Beq)
    kst = keystoneness(m, Beq)
    @test kst == abs.([S[2, 1], S[1, 2]])
    # Three species with competition.
    fw = Foodweb(zeros(Bool, 3, 3))
    m = default_model(fw, ProducersCompetition([1 0.1 0.2; 0.2 1 0.1; 0.03 0.02 1]))
    Beq = simulate(m, [1, 1, 1], 1_000; show_degenerated).u[end]
    S = sensitivity(m, Beq)
    kst = keystoneness(m, Beq)
    @test kst == [
        abs(S[2, 1]) + abs(S[3, 1]),
        abs(S[1, 2]) + abs(S[3, 2]),
        abs(S[1, 3]) + abs(S[2, 3]),
    ]
end

@testset "Resistance." begin
    # Two independent species.
    fw = Foodweb([0 0; 0 0])
    m = default_model(fw)
    Beq = simulate(m, [1, 1], 100; show_degenerated = false).u[end]
    @test resistance(m, Beq) == [-1, -1]
    @test resistance(m, Beq; perturbation_on = 1, response_of = 2) == 0
    @test resistance(m, Beq; perturbation_on = 2, response_of = 1) == 0
    @test resistance(m, Beq; perturbation_on = 1, response_of = 1) == -1
    @test resistance(m, Beq; perturbation_on = 2, response_of = 2) == -1
    @test resistance(m, Beq; aggregated = true) == -2
    @test resistance(m, Beq; perturbation_on = 1, aggregated = true) == -1
    @test resistance(m, Beq; response_of = 1, aggregated = true) == -1
    # Consumer - resource.
    fw = Foodweb([2 => 1])
    m = default_model(fw)
    Beq = simulate(m, [1, 1], 100; show_degenerated = false).u[end]
    @test isapprox(resistance(m, Beq; perturbation_on = 1, response_of = 1), 0, atol = 1e-6)
    @test resistance(m, Beq; perturbation_on = 1, response_of = 2) < 0
    @test resistance(m, Beq; perturbation_on = 2, response_of = 1) > 0
    # Test keyword arguments.
    res_exp = resistance(m, Beq)
    @test resistance(m, Beq; perturbation_on = [1, 2]) == res_exp
    @test resistance(m, Beq; perturbation_on = [2, 1]) == res_exp
    @test resistance(m, Beq; response_of = [1, 2]) == res_exp
    @test resistance(m, Beq; response_of = [2, 1]) == res_exp[[2, 1]]
end

@testset "Resistance from simulations." begin
    # Two independent species.
    fw = Foodweb([0 0; 0 0])
    m = default_model(fw)
    Beq = simulate(m, [1, 1], 100; show_degenerated = false).u[end]
    res = resistance_simulation(m, Beq)
    @test isapprox(res, [-1, -1]; atol = 1e-6)
    res = resistance_simulation(m, Beq; mortality_increment = [0.1, 0])
    @test isapprox(res[1], -1, atol = 1e-2)
    @test res[2] == :undefined
    res = resistance_simulation(m, Beq; mortality_increment = [0.1, 0], normalized = false)
    @test res[2] == 0
    @test isapprox(resistance_simulation(m, Beq; aggregated = true), -2, atol = 1e-2)
end
