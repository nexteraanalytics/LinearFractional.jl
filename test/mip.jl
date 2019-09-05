@testset "Binary Variables" begin
    lfp = LinearFractionalModel(with_optimizer(Cbc.Optimizer))

    x1 = @variable(lfp, lower_bound=0, base_name="x1", binary=true)
    x2 = @variable(lfp, lower_bound=0, binary=true, base_name="x2")
    @constraint(lfp, -x1 + x2 <= 4)
    @constraint(lfp, 2x1 + x2 <= 14)

    LinearFractional.set_objective(lfp, MOI.MIN_SENSE,
              @expression(lfp, -2x1 + x2 + 2),
              @expression(lfp, x1 + -x2 + 4))

    optimize!(lfp)

    @test termination_status(lfp) === MOI.OPTIMAL
    @test value(x1) ≈ 1.0
    @test value(x2) ≈ 0.0

    @testset "Set the binary big-M" begin
        lfp = LinearFractionalModel(with_optimizer(Cbc.Optimizer); binary_M = 100.0)
        @test lfp.options.binary_M == 100.0

        lfp = LinearFractionalModel(with_optimizer(Cbc.Optimizer); binary_M = 10.0)
        @test lfp.options.binary_M == 10.0
    end
end
