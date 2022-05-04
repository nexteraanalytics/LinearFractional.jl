# Approach based on JuMP/test/JuMPExtension.jl
module LinearFractional

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
using JuMP
import JuMP: @constraint,
    @variable,
    AbstractModel,
    AbstractVariableRef,
    AbstractJuMPScalar,
    AffExpr,
    Model,
    ConstraintRef,
    GenericAffExpr,
    GenericQuadExpr,
    VariableRef,
    VariableInfo,
    _VariableInfoExpr,
    add_to_expression!,
    objective_function,
    optimize!,
    set_objective,
    lower_bound,
    has_lower_bound,
    set_lower_bound,
    delete_lower_bound,
    has_upper_bound,
    set_upper_bound,
    delete_upper_bound,
    upper_bound,
    is_fixed,
    fix_value,
    fix,
    unfix,
    start_value,
    set_start_value,
    set_binary,
    unset_binary,
    is_binary,
    is_integer,
    set_integer,
    unset_integer,
    name,
    set_name,
    owner_model,
    linear_terms,
    index,
    Containers,
    ScalarVariable,
    add_variable,
    constant,
    termination_status,
    objective_value,
    result_count

import Base.convert

using MacroTools

export LinearFractionalModel

macro forward(ex, fs)
    T, field = ex.args[1], ex.args[2].value

    T = esc(T)
    fs = isexpr(fs, :tuple) ? map(esc, fs.args) : [esc(fs)]
    :($([:($f(x::$T, args...) = (Base.@_inline_meta; $f(x.$field, args...)))
        for f in fs]...);
    nothing)
end

Base.@kwdef struct LinearFractionalOptions
    binary_M::Float64 = 1e3
end

mutable struct LinearFractionalModel <: AbstractModel
    model::JuMP.Model #JuMP.Model(solver=solver)
    t::VariableRef
    denominator_constraint::Union{Nothing, ConstraintRef}
    options::LinearFractionalOptions

    function LinearFractionalModel(; kwargs...)
        all_kwargs = collect(kwargs)
        jump_kwargs = filter(pair -> pair.first ∉ fieldnames(LinearFractionalOptions),
            all_kwargs)
        lf_option_kwargs = filter(pair -> pair.first ∈ fieldnames(LinearFractionalOptions),
            all_kwargs)

        model = JuMP.Model(jump_kwargs...)
        t = @variable(model, lower_bound = 10*eps(), base_name = "lf_aux")

        options = LinearFractionalOptions(; lf_option_kwargs...)
        new(model, t, nothing, options)
    end
end

function LinearFractionalModel(optimizer;
               add_bridges::Bool=true, kwargs...)
    model = LinearFractionalModel(; kwargs...)
    JuMP.set_optimizer(model.model, optimizer,
                  add_bridges=add_bridges)
    return model
end

@forward LinearFractionalModel.model JuMP.constraint_by_name,
    JuMP.delete,
    JuMP.is_valid,
    JuMP.num_variables,
    JuMP.objective_function,
    JuMP.objective_function_type,
    JuMP.objective_sense,
    JuMP.objective_value,
    JuMP.optimize!,
    JuMP.set_objective_sense,
    JuMP.termination_status,
    JuMP.variable_by_name,
    JuMP.result_count

Base.broadcastable(model::LinearFractionalModel) = Ref(model)

JuMP.object_dictionary(model::LinearFractionalModel) = model.model.obj_dict

# Wrap `VariableRef` instead of having fields, say
# model::LinearFractionalModel and idx::Int
# so that we can leverage existing `add_variable` and other methods
# by simply passing the `vref`
struct LinearFractionalVariableRef <: JuMP.AbstractVariableRef
    model::LinearFractionalModel
    vref::VariableRef
end

@forward LinearFractionalVariableRef.vref name,
    set_name,
    lower_bound,
    has_lower_bound,
    set_lower_bound,
    delete_lower_bound,
    has_upper_bound,
    set_upper_bound,
    delete_upper_bound,
    upper_bound,
    is_fixed,
    fixed_value,
    fix,
    unfix,
    start_value,
    set_start_value,
    is_binary,
    is_integer,
    set_integer,
    unset_integer,
    owner_model,
    _assert_isfinite

Base.copy(v::LinearFractionalVariableRef) = v
Base.copy(v::LinearFractionalVariableRef, new_model::LinearFractionalModel) = LinearFractionalVariableRef(new_model, v.idx)

Base.:(==)(v::LinearFractionalVariableRef, w::LinearFractionalVariableRef) = v.model === w.model && (v.vref == w.vref)
Base.broadcastable(v::LinearFractionalVariableRef) = Ref(v)
JuMP.isequal_canonical(v::LinearFractionalVariableRef, w::LinearFractionalVariableRef) = v == w
function JuMP.add_variable(m::LinearFractionalModel, v::JuMP.AbstractVariable, name::String="")

    inner_vref = JuMP.add_variable(m.model, v, name)
    vref = LinearFractionalVariableRef(m, inner_vref)

    if has_lower_bound(inner_vref)
        lb = lower_bound(inner_vref)
        delete_lower_bound(inner_vref)
        set_lower_bound(vref, lb)
    end

    if has_upper_bound(inner_vref)
        ub = upper_bound(inner_vref)
        delete_upper_bound(inner_vref)
        set_upper_bound(vref, ub)
    end

    if is_binary(inner_vref)
        unset_binary(inner_vref)
        set_binary(vref)
    end

    return vref

end
function JuMP.delete(model::LinearFractionalModel, vref::LinearFractionalVariableRef)
    JuMP.delete(model.model, vref.vref)
end
function JuMP.is_valid(model::LinearFractionalModel, vref::LinearFractionalVariableRef)
    return (model === vref.model &&
            is_valid(model.model, vref.vref))
end


# Internal functions

# Defined just to be forwarded above
variable_info(vref::VariableRef) = vref.model.variables[vref.idx].info
function update_variable_info(vref::VariableRef, info::JuMP.VariableInfo)
   vref.model.variables[vref.idx] = JuMP.ScalarVariable(info)
end

function set_lower_bound(vref::LinearFractionalVariableRef, lower::Number)
    @constraint(vref.model, vref >= lower)
end

function set_upper_bound(vref::LinearFractionalVariableRef, upper::Number)
    @constraint(vref.model, vref <= upper)
end

"""
If the original variable is binary,
then y ∈ {0, t}.  We can solve this with a big-M formuation:
e.g., y = zt for a binary variable z ∈ {0, 1}

y ∈ [0, t]
y >= 0   (1)
y <= t   (2)

if z = 0, y <= 0 --> y == 0
if z = 1, y is free
y <= zM  (3)

if z = 0, y free
if z = 1, y >= t --> y == t
y >= t - (1-z)M   (4)
"""
function JuMP.set_binary(lvref::LinearFractionalVariableRef)
    M = lvref.model.options.binary_M

    trans_model = lvref.model.model
    t = lvref.model.t
    y = lvref.vref
    z = @variable(trans_model, binary=true)

    # 0 ≦ y ≦ t
    set_lower_bound(y, 0.0)
    @constraint(trans_model, y <= t)

    # y <= zM  (3)
    @constraint(trans_model, y <= z * M)

    # y >= t - (1-z)M   (4)
    @constraint(trans_model, y >= t - (1 - z) * M)
end


function transform_constraint(model::LinearFractionalModel, constraint_ref::ConstraintRef)
    α = convert(Float64, normalized_rhs(constraint_ref))
    set_normalized_rhs(constraint_ref, 0.0)
    JuMP.set_normalized_coefficient(constraint_ref, model.t, -α)
end

function JuMP.add_constraint(model::LinearFractionalModel, c::JuMP.AbstractConstraint,
                             name::String="")
    cref = JuMP.add_constraint(model.model, c, name)
    transform_constraint(model, cref)
    return cref
end

# Objective
function transform_objective_numerator(m::LinearFractionalModel)
    obj = objective_function(m.model, AffExpr)
    β = obj.constant
    obj.constant = 0.0
    JuMP.add_to_expression!(obj, β, m.t)
    JuMP.set_objective_function(m.model, obj)
end

function set_numerator(m::LinearFractionalModel, sense::MOI.OptimizationSense, f)
    JuMP.set_objective(m.model, sense, f)
    transform_objective_numerator(m)
end

function JuMP.set_objective(m::LinearFractionalModel, sense::MOI.OptimizationSense,
    numerator::GenericAffExpr{Float64, LinearFractionalVariableRef},
    denominator::GenericAffExpr{Float64, LinearFractionalVariableRef})
    set_numerator(m, sense, numerator)
    set_denominator(m, denominator)
end

function JuMP.set_objective(m::LinearFractionalModel, sense::MOI.OptimizationSense,
    numerator::GenericAffExpr{Float64, LinearFractionalVariableRef},
    denominator::Float64)

    denominator <= 0.0 && error("A LinearFractionalModel with a constant denominator must have a positive denominator.")

    set_numerator(m, sense, numerator)
    @constraint(m.model, m.t == 1.0/denominator)
end



"""
Given a denominator

∑d_j x_i + β,

the transformed problem has the constraint

∑d_j x_i + βt == 1.

This cannot be achieved by adding a pseudo-constraint to the original problem
because of the existence of the constant in the transformed constraint.
Normal transformed constraints have 0 on the RHS and no constant.  Therefore,
we need to convert the GenericAffExpr{Float64,LinearFractionalVariableRef} to
the corresponding GenericAffExpr{Float64, VariableRef} in order to create the
final expression in the denominator.
"""
function set_denominator(model::LinearFractionalModel, expr_lf)
    if !isnothing(model.denominator_constraint)
        delete(model, model.denominator_constraint)
        model.denominator_constraint = nothing
    end
    β = expr_lf.constant
    expr = convert(GenericAffExpr{Float64,VariableRef}, expr_lf)
    expr.constant = 0.0
    add_to_expression!(expr, β, model.t)
    model.denominator_constraint = @constraint(model.model, expr == 1.0)
end

"""
Convert a GenericAffExpr of LinearFractionalVariableRef to a GenericAffExpr of VariableRef.
Note that this does not perform the constant -> t x constant transformation since this
is taken care of in add_constaint.  The separation allows this transformation to
be used in `set_denominator`
"""
function convert(::Type{GenericAffExpr{Float64,VariableRef}}, expr_lf::GenericAffExpr{Float64,LinearFractionalVariableRef})
    GenericAffExpr{Float64,VariableRef}(expr_lf.constant, [Pair(term.vref, coef) for (term, coef) in expr_lf.terms])
end




# Show
function JuMP.show_backend_summary(io::IO, model::LinearFractionalModel) end
function JuMP.show_objective_function_summary(io::IO, model::LinearFractionalModel)
    println(io, "Objective function type: ",
            JuMP.objective_function_type(model.model))
end
function JuMP.objective_function_string(print_mode, model::LinearFractionalModel)
    return JuMP.function_string(print_mode, JuMP.objective_function(model.model))
end
_plural(n) = (isone(n) ? "" : "s")
function JuMP.show_constraints_summary(io::IO, model::LinearFractionalModel)
    JuMP.show_constraints_summary(io, model.model)
end
function JuMP.constraints_string(print_mode, model::LinearFractionalModel)
    JuMP.constraints_string(print_mode, model.model)
end

function JuMP.value(v::LinearFractionalVariableRef; result=1)::Float64
    return JuMP.value(v.vref, result=result)/JuMP.value(v.model.t, result=result)
end


# These functions that feel like they break the LinearFractionalVariableRef
# API.
function _assert_isfinite(a::GenericAffExpr{Float64,LinearFractionalVariableRef})
    for (coef, var) in linear_terms(a)
        isfinite(coef) || error("Invalid coefficient $coef on variable $var.")
    end
end


function MOI.ScalarAffineFunction(a::GenericAffExpr{Float64,LinearFractionalVariableRef})
    _assert_isfinite(a)
    terms = MOI.ScalarAffineTerm{Float64}[MOI.ScalarAffineTerm(t[1],
                                                               index(t[2].vref))
                                          for t in linear_terms(a)]
    return MOI.ScalarAffineFunction(terms, a.constant)
end


end
