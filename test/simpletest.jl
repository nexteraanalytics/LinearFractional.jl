# Validated from http://www.ams.jhu.edu/~castello/625.414/Handouts/FractionalProg.pdf

@testset "Simple example" begin
    lfp = LinearFractionalModel(with_optimizer(Clp.Optimizer))
    x1 = @variable(lfp, lower_bound=0, base_name="x1")
    x2 = @variable(lfp, lower_bound=0, upper_bound=6, base_name="x2")
    @constraint(lfp, -x1 + x2 <= 4)
    @constraint(lfp, 2x1 + x2 <= 14)

    LinearFractional.set_objective(lfp, MOI.MIN_SENSE,
              @expression(lfp, -2x1 + x2 + 2),
              @expression(lfp, x1 + 3x2 + 4))

    optimize!(lfp)
    @test termination_status(lfp) === MOI.OPTIMAL
    @test value(x1) ≈ 7.0
    @test value(x2) ≈ 0.0
end
