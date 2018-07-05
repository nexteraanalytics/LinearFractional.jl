__precompile__()

module LinearFractional

using JuMP
using JuMP
import JuMP:
    getobjectivevalue,
    getvalue,
    addconstraint,
    solve,
    validmodel,
    constructvariable!,
    constructconstraint!,
    variabletype,
    JuMPArray,
    JuMPContainerData,
    setlowerbound,
    setupperbound,
    addtoexpr_reorder,
    assert_validmodel,
    setobjective,
    storecontainerdata
import MathProgBase.AbstractMathProgSolver
importall Base.Operators
using Parameters

export LinearFractionalModel,
    @denominator,
    @numerator

@with_kw mutable struct LinearFractionalModel <: JuMP.AbstractModel
    transformedmodel::JuMP.Model
    t::JuMP.Variable
    denom::AffExpr
    dictList::Vector=Any[]
end


function LinearFractionalModel(;solver::AbstractMathProgSolver=ClpSolver())
    model = JuMP.Model(solver=solver)
    t = @variable(model, lowerbound=eps(Float64), basename="t")
    denom = AffExpr()
    LinearFractionalModel(model, t, denom, Any[])
end


struct LinearFractionalVariable <: JuMP.AbstractJuMPScalar
    ## Variable in the untransformed space
    model::LinearFractionalModel
    var::JuMP.Variable ## Internal variable in the transformed space
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

# Need to define Base.one to get dot() for free, but I don't think we can
# do this without knowing the model (to get t)
# Base.one(::Type{LinearFractionalVariable}) = LinearFractionalAffExpr(Variable[],Float64[],1.0)


storecontainerdata(m::LinearFractionalModel, variable, varname, idxsets, idxpairs, condition) =
    m.transformedmodel.varData[variable] = JuMPContainerData(varname, map(collect,idxsets), idxpairs, condition)

struct LinearFractionalAffExpr
    afftrans::AffExpr
    t::JuMP.Variable
end



function AffExpr(vars::Vector{LinearFractionalVariable}, coeffs, constant)
    t = vars[1].model.t
    ys = [var.var for var in vars]
    return LinearFractionalAffExpr(AffExpr(ys, coeffs, 0) + constant * t, t)
end




function addvariable(model::LinearFractionalModel, lb::Number, ub::Number, basename::String)
    var = @variable(model.transformedmodel, basename=basename)
    if !isinf(lb) || !isinf(ub)
        addconstraint(model, lb, ub)
    end
    return LinearFractionalVariable(model, var)
end


addvariable(model::LinearFractionalModel, basename::String) = addvariable(model, -Inf, Inf, basename)





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
include("constraints.jl")
include("macros.jl")
include("parseexpr.jl")

end
