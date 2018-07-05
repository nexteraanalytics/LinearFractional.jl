@testset "Array variables" begin
    lfp = LinearFractionalModel(solver=ClpSolver())
    a = [-2, 0]
    x = @variable(lfp, [i=1:2], basename="x", lowerbound=0)
    @constraint(lfp, x[2] <= 6)
    @constraint(lfp, -x[1] + x[2] <= 4)
    @constraint(lfp, 2x[1] + x[2] <= 14)
    @numerator(lfp,  :Min, sum(a[i] * x[i] for i in 1:2) + 2)
    @denominator(lfp,  x[1] + 3x[2] + 4)
    @test solve(lfp) == :Optimal
    @test all(getvalue(x) .≈ [7.0, 0.0])
end


@testset "Array constraints" begin
    lfp = LinearFractionalModel(solver=ClpSolver())
    x = @variable(lfp, [i=1:2], basename="x")
    a = [2, 1]
    upbs = [4, 20]
    lbs = [-4, -20]
    @constraint(lfp, [i=1:2], x[i] <= upbs[i])
    @constraint(lfp, [i=1:2], x[i] >= lbs[i])
    @numerator(lfp,  :Min, sum(a[i] * x[i] for i in 1:2) + 2)
    @denominator(lfp,  -x[1] + 20*x[2] + 4)
    @test solve(lfp) == :Optimal
    @test all(getvalue(x) .≈ [-4.0, -0.2])

    # (sum(a[i] * getvalue(x)[i] for i in 1:2) + 2) / (-getvalue(x)[1] + 20*getvalue(x)[2] + 4)
    # getobjectivevalue(lfp)
    #
    # f(x) = (sum(a[i] * x[i] for i in 1:2) + 2) / (-x[1] + 20*x[2] + 4)
end
