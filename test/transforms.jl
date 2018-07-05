@testset "Transformations" begin
    lfm = LinearFractionalModel(solver=ClpSolver())
    x = @variable(lfm)
    xplus5 = x + 5.0
    @test xplus5.afftrans.vars == [x.var, lfm.t]
    @test xplus5.afftrans.coeffs == [1.0, 5.0]
    @test xplus5.afftrans.constant == 0.0

    xtimes2 = 2.0 * x
    @test xtimes2.afftrans.vars == [x.var]
    @test xtimes2.afftrans.coeffs == [2.0]
    @test xtimes2.afftrans.constant == 0.0
end
