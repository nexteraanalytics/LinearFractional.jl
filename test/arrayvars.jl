@testset "Array variables" begin
    lfp = LinearFractionalModel(with_optimizer(Clp.Optimizer))
    a = [-2, 0]
    x = @variable(lfp, [i=1:2], base_name="x", lower_bound=0)
    @constraint(lfp, x[2] <= 6)
    @constraint(lfp, -x[1] + x[2] <= 4)
    @constraint(lfp, 2x[1] + x[2] <= 14)
    numer = @expression(lfp, sum(a[i] * x[i] for i in 1:2) + 2)
    denom = @expression(lfp, x[1] + 3x[2] + 4)
    set_objective(lfp, MOI.MIN_SENSE, numer, denom)
    optimize!(lfp)
    @test termination_status(lfp) == MOI.OPTIMAL
    @test all(value.(x) .≈ [7.0, 0.0])
end


@testset "Array constraints" begin
    lfp = LinearFractionalModel(with_optimizer(Clp.Optimizer))
    x = @variable(lfp, [i=1:2], base_name="x")
    a = [2.0, 1.0]
    upbs = [4.0, 20.0]
    lbs = [-4.0, -20.0]
    @constraint(lfp, [i=1:2], x[i] <= upbs[i])
    @constraint(lfp, [i=1:2], x[i] >= lbs[i])
    numer = @expression(lfp,  sum(a[i] * x[i] for i in 1:2) + 2)
    denom = @expression(lfp,  -x[1] + 20*x[2] + 4)
    set_objective(lfp, MOI.MIN_SENSE, numer, denom)
    optimize!(lfp)
    @test termination_status(lfp) === MOI.OPTIMAL
    @test value.(x) ≈ [-4.0, -0.2]
end


@testset "Match LP" begin
    lfp = LinearFractionalModel(with_optimizer(Clp.Optimizer))
    a = [4, 2]
    upbs = [4.0, 20.0]
    lbs = [-1.0, -10.0]
    x = @variable(lfp, [i=1:2], base_name="x", lower_bound=lbs[i], upper_bound=upbs[i])
    @constraint(lfp, x[1] + x[2] <= 5.0)
    @constraint(lfp, x[1] - 2*x[2] >= 10.0)
    numer = @expression(lfp,  dot(a, x))
    denom = @expression(lfp,  sum(x))
    set_objective(lfp, MOI.MIN_SENSE, numer, denom)
    optimize!(lfp)
    m = Model(with_optimizer(Clp.Optimizer))
    xm = @variable(m, [i=1:2], base_name="x")
    t = @variable(m, lower_bound=0.0)
    a = [4, 2]
    upbs = [4, 20]
    lbs = [-1,-10]
    @constraint(m, xm[1] + xm[2] <= 5.0t)
    @constraint(m, xm[1] - 2*xm[2] >= 10.0t)
    @constraint(m, [i=1:2], xm[i] <= upbs[i] * t)
    @constraint(m, [i=1:2], xm[i] >= lbs[i] * t)
    @constraint(m, sum(xm) == 1)
    @objective(m,  Min, sum(a[i] * xm[i] for i in 1:2))
    optimize!(m)
    @test value.(xm)/value(t) ≈ value.(x)
    @test objective_value(lfp) ≈ objective_value(m)
end



@testset "Constant denominator" begin
    lfp = LinearFractionalModel(with_optimizer(Clp.Optimizer))
    a = [4, 2]
    upbs = [4, 20]
    lbs = [-1, -10]
    x = @variable(lfp, [i=1:2], base_name="x", lower_bound=lbs[i], upper_bound=upbs[i])
    @constraint(lfp, x[1] + x[2] <= 5.0)
    @constraint(lfp, x[1] - 2*x[2] >= 10.0)
    @constraint(lfp, [i=1:2], x[i] <= upbs[i])
    @constraint(lfp, [i=1:2], x[i] >= lbs[i])
    numer = @expression(lfp, sum(a[i] * x[i] for i in 1:2))
    denom = @expression(lfp, 5.0)
    set_objective(lfp, MOI.MIN_SENSE, numer, denom)
    optimize!(lfp)

    m = Model(with_optimizer(Clp.Optimizer))
    a = [4, 2]
    upbs = [4, 20]
    lbs = [-1,-10]
    xm = @variable(m, [i=1:2], base_name="x", lower_bound=lbs[i], upper_bound=upbs[i])
    @constraint(m, xm[1] + xm[2] <= 5.0)
    @constraint(m, xm[1] - 2*xm[2] >= 10.0)
    @objective(m,  Min, sum(a[i] * xm[i] for i in 1:2)/5.0)
    optimize!(m)
    @test value.(xm) ≈ value.(x)
    @test objective_value(lfp) ≈ objective_value(m)
end
