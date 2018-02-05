module LinearFractional

using JuMP
using JuMP
import JuMP: getobjectivevalue, getvalue, addconstraint, solve, validmodel,
    constructvariable!, constructconstraint!, variabletype, JuMPArray,
    setlowerbound, setupperbound, addtoexpr_reorder, assert_validmodel,
    setobjective
importall Base.Operators
using Parameters

export LinearFractionalModel,
    @denominator,
    @numerator

mutable struct LinearFractionalModel
    transformedmodel::JuMP.Model
    t::JuMP.Variable
    denom::AffExpr
end


function LinearFractionalModel(;solver=ClpSolver())
    model = JuMP.Model(solver=solver)
    t = @variable(model, lowerbound=eps(Float64), basename="t")
    denom = AffExpr()
    LinearFractionalModel(model, t, denom)
end


struct LinearFractionalVariable
    model::LinearFractionalModel
    var::JuMP.Variable
end


function LinearFractionalVariable(m::LinearFractionalModel, lower::Number, upper::Number, cat::Symbol, name::AbstractString="", value::Number=NaN)
    var = JuMP.Variable(m.transformedmodel, -Inf, Inf, cat, name, value)
    lfvar = LinearFractionalVariable(m, var)
    if !isinf(lower)
        setlowerbound(lfvar, lower)
    end
    if !isinf(upper)
        setupperbound(lfvar, upper)
    end
    return lfvar
end

LinearFractionalVariable(m::Model,lower::Number,upper::Number,cat::Symbol,objcoef::Number,
    constraints::JuMPArray,coefficients::Vector{Float64}, name::AbstractString="", value::Number=NaN) =
    LinearFractionalVariable(m, lower, upper, cat, objcoef, constraints.innerArray, coefficients, name, value)

getvalue(x::LinearFractionalVariable) = getvalue(x.var)/getvalue(x.model.t)


struct LinearFractionalAffExpr
    afftrans::AffExpr
    t::JuMP.Variable
end


struct LinearFractionalLinearConstraint
    lctrans::Vector{LinearConstraint}  # 1 or 2 (lb, ub)
end


function AffExpr(vars::Vector{LinearFractionalVariable}, coeffs, constant)
    t = vars[1].model.t
    ys = [var.var for var in vars]
    return LinearFractionalAffExpr(AffExpr(ys, coeffs, 0) + constant * t, t)
end


function LinearConstraint(aff::LinearFractionalAffExpr, lb, ub)
    cons = Vector{LinearConstraint}()
    if lb == ub
         push!(cons, LinearConstraint(aff.afftrans - aff.t * lb, 0, 0))
    else
        if !isinf(lb)
            push!(cons, LinearConstraint(aff.afftrans - aff.t * lb, 0, Inf))
        end
        if !isinf(ub)
            push!(cons, LinearConstraint(aff.afftrans - aff.t * ub, -Inf, 0))
        end
    end
    LinearFractionalLinearConstraint(cons)
end


function addvariable(model::LinearFractionalModel, lb::Number, ub::Number, basename::String)
    var = @variable(model.transformedmodel, basename=basename)
    if !isinf(lb) || !isinf(ub)
        addconstraint(model, lb, ub)
    end
    return LinearFractionalVariable(model, var)
end


addvariable(model::LinearFractionalModel, basename::String) = addvariable(model, -Inf, Inf, basename)


function addconstraint(model::LinearFractionalModel, constraint::LinearFractionalLinearConstraint)
    for con in constraint.lctrans
        addconstraint(model.transformedmodel, con)
    end
end


function setdenominator!(m::LinearFractionalModel, aff::LinearFractionalAffExpr)
    addconstraint(m.transformedmodel, LinearConstraint(aff.afftrans, 1, 1))
    m.denom = aff.afftrans
end


function setobjective(m::LinearFractionalModel, sense::Symbol, numer::LinearFractionalAffExpr, denom::LinearFractionalAffExpr)
    JuMP.setobjective(m.transformedmodel, sense, numer.afftrans)
    setdenominator!(m, denom)
end


function solve(model::LinearFractionalModel)
    JuMP.solve(model.transformedmodel)
end


getobjectivevalue(model::LinearFractionalModel) = getobjectivevalue(model.transformedmodel)/getvalue(model.denom)


validmodel(model::LinearFractionalModel, name) = nothing


function setlowerbound(v::LinearFractionalVariable, lb)
    if iszero(lb)
        setlowerbound(v.var, 0.0)
    else
        addconstraint(v.model, LinearConstraint(AffExpr([v], [1], 0), lb, Inf))
    end
end


function setupperbound(v::LinearFractionalVariable, ub)
    if iszero(ub)
        setupperbound(v.var, 0.0)
    else
        addconstraint(v.model, LinearConstraint(AffExpr([v], [1], 0), -Inf, ub))
    end
end


setname(v::LinearFractionalVariable, name) = setname(v.var, name)

function setobjective(m::LinearFractionalModel, sense::Symbol, numer::LinearFractionalAffExpr, denom::LinearFractionalAffExpr)
    if sense == :Min
        moisense = MOI.MinSense
    else
        @assert sense == :Max
        moisense = MOI.MaxSense
    end
    MOI.set!(m.transformedmodel.moibackend, MOI.ObjectiveSense(), moisense)
    MOI.set!(m.transformedmodel.moibackend, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction(numer.afftrans))
    setdenominator!(m, denom)
    nothing
end

include("operators.jl")
include("macros.jl")
include("parseexpr.jl")

end
