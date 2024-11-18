@testset "Basic SOCP with linear fractional objective" begin
    model = LinearFractionalModel(ECOS.Optimizer)
    
    @variable(model, y₁ >= 1)
    @variable(model, y₂)
    @constraint(model, y₁ <= 2)
    @constraint(model, [1y₁; y₂] in SecondOrderCone())
    
    set_objective(model, MOI.MIN_SENSE, y₁ + 2y₂, y₁ + 1)
    optimize!(model)
    
    @test termination_status(model) === MOI.OPTIMAL

    # worked out manually
    @test isapprox(value(y₁), 2.0, rtol=1e-4)
    @test isapprox(value(y₂), -2.0, rtol=1e-4)
    @test isapprox(objective_value(model), -2/3, rtol=1e-4)
end


@testset "Multiple cone constraints with optimality verification" begin
    # Solve the original problem
    model = LinearFractionalModel(ECOS.Optimizer)
    
    @variable(model, x[1:2] >= 0.1)
    
    set_objective(model, MOI.MIN_SENSE, 
        2x[1] + x[2] + 1,   # numerator
        x[1] + 2x[2] + 2    # denominator
    )

    @constraint(model, [1.0; x[1]; x[2]] in SecondOrderCone())
    @constraint(model, [x[1] + x[2]; x[1]; x[2]] in RotatedSecondOrderCone())

    optimize!(model)
    
    @test termination_status(model) === MOI.OPTIMAL
    opt_x = value.(x)
    opt_val = objective_value(model)
    
    # Test cone constraints satisfaction
    @test norm([opt_x[1]; opt_x[2]]) <= 1.0 + 1e-4
    @test 2 * (opt_x[1] + opt_x[2]) * opt_x[1] >= opt_x[2]^2 - 1e-4
    @test all(opt_x .>= 0.1 - 1e-4)

    # Verify optimality by checking feasibility of perturbed problems

    function solve_feasibility_problem(γ)
        m = Model(ECOS.Optimizer)
        
        @variable(m, x[1:2] >= 0.1)
        @constraint(m, [1.0; x[1]; x[2]] in SecondOrderCone())
        @constraint(m, [x[1] + x[2]; x[1]; x[2]] in RotatedSecondOrderCone())
        
        @constraint(m, 2x[1] + x[2] + 1 <= γ * (x[1] + 2x[2] + 2))
        
        optimize!(m)
        return termination_status(m)
    end

    # Test with perturbations around optimal value
    @test solve_feasibility_problem(opt_val + 1e-4) == MOI.OPTIMAL  # Should be feasible
    @test solve_feasibility_problem(opt_val - 1e-4) == MOI.INFEASIBLE  # Should be infeasible
end
