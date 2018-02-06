@testset "Try a problem with array variables" begin
    lfp = LinearFractionalModel(solver=ClpSolver())
    x = @variable(lfp, [i=1:2], basename="x", lowerbound=0)
    @constraint(lfp, x[2] <= 6)
    @constraint(lfp, -x[1] + x[2] <= 4)
    @constraint(lfp, 2x[1] + x[2] <= 14)
    @numerator(lfp,  :Min, -2x[1] + x[2] + 2)
    @denominator(lfp,  x[1] + 3x[2] + 4)
    @test solve(lfp) == :Optimal
    @test all(getvalue(x) .≈ [7.0, 0.0])
end
