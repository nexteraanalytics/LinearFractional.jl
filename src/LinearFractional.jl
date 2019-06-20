# Approach based on JuMP/test/JuMPExtension.jl
module LinearFractional

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
import JuMP
import JuMP: @constraint,
    @variable,
    AbstractVariableRef,
    AffExpr,
    ConstraintRef,
    GenericAffExpr,
    GenericQuadExpr,
    VariableRef,
    add_to_expression!,
    objective_function,
    optimize!,
    set_objective

using MacroTools

export @denominator,
    @numerator,
    LinearFractionalModel

include("parseexpr.jl")

macro forward(ex, fs)
    T, field = ex.args[1], ex.args[2].value

    T = esc(T)
    fs = isexpr(fs, :tuple) ? map(esc, fs.args) : [esc(fs)]
    :($([:($f(x::$T, args...) = (Base.@_inline_meta; $f(x.$field, args...)))
        for f in fs]...);
    nothing)
end

# struct ConstraintIndex
#     value::Int # Index in `model.constraints`
# end

mutable struct LinearFractionalModel <: JuMP.AbstractModel
    model::JuMP.Model #JuMP.Model(solver=solver)
    t::VariableRef
    denominator_constraint::Union{Nothing, ConstraintRef}
    #
    # nextvaridx::Int                                 # Next variable index is nextvaridx+1
    # variables::Dict{Int, JuMP.ScalarVariable}       # Map varidx -> variable
    # var_to_name::Dict{Int, String}                  # Map varidx -> name
    # name_to_var::Union{Dict{String, Int}, Nothing}  # Map varidx -> name
    # nextconidx::Int                                 # Next constraint index is nextconidx+1
    # constraints::Dict{ConstraintIndex,
    #                   JuMP.AbstractConstraint}      # Map conidx -> variable
    # con_to_name::Dict{ConstraintIndex, String}      # Map conidx -> name
    # name_to_con::Union{Dict{String, ConstraintIndex},
    #                    Nothing}                     # Map name -> conidx
    # objectivesense::MOI.OptimizationSense
    # objective_function::JuMP.AbstractJuMPScalar
    # obj_dict::Dict{Symbol, Any}                     # Same that JuMP.Model's field `obj_dict`
    function LinearFractionalModel(; kwargs...)
        model = JuMP.Model(kwargs...)
        t = @variable(model, lower_bound = 0.0, base_name = "lf_aux")
        new(model, t, nothing)
        # new(0, Dict{Int, JuMP.AbstractVariable}(),
        #     Dict{Int, String}(), nothing,                        # Variables
        #     0, Dict{ConstraintIndex, JuMP.AbstractConstraint}(),
        #     Dict{ConstraintIndex, String}(), nothing,            # Constraints
        #     MOI.FEASIBILITY_SENSE,
        #     zero(JuMP.GenericAffExpr{Float64, MyVariableRef}),
        #     Dict{Symbol, Any}())
    end
end

function LinearFractionalModel(optimizer_factory::JuMP.OptimizerFactory;
               bridge_constraints::Bool=true, kwargs...)
    model = LinearFractionalModel(; kwargs...)
    JuMP.set_optimizer(model.model, optimizer_factory,
                  bridge_constraints=bridge_constraints)
    return model
end

Base.broadcastable(model::LinearFractionalModel) = Ref(model)

JuMP.object_dictionary(model::LinearFractionalModel) = model.model.obj_dict

# TODO:  Maybe we do need to keep our wrapper variable type??  Otherwise
# value(x) will be the raw (not back-transformed version)
# also, and more importantly, variable bounds are not transformed...


# # Variables
# struct MyVariableRef <: JuMP.AbstractVariableRef
#     model::LinearFractionalModel # `model` owning the variable
#     idx::Int       # Index in `model.variables`
# end
# Base.copy(v::MyVariableRef) = v
Base.copy(v::VariableRef, new_model::LinearFractionalModel) = VariableRef(new_model, v.idx)
#
# Base.:(==)(v::MyVariableRef, w::MyVariableRef) = v.model === w.model && v.idx == w.idx
# Base.broadcastable(v::MyVariableRef) = Ref(v)
# JuMP.isequal_canonical(v::MyVariableRef, w::MyVariableRef) = v == w
JuMP.variable_type(::LinearFractionalModel) = VariableRef
function JuMP.add_variable(m::LinearFractionalModel, v::JuMP.AbstractVariable, name::String="")
    JuMP.add_variable(m.model, v, name)
    # m.nextvaridx += 1
    # vref = MyVariableRef(m, m.nextvaridx)
    # m.variables[vref.idx] = v
    # JuMP.set_name(vref, name)
    # vref
end
function JuMP.delete(model::LinearFractionalModel, vref::VariableRef)
    @assert JuMP.is_valid(model.model, vref)
    delete!(model.model.variables, vref.idx)
    delete!(model.model.var_to_name, vref.idx)
end
function JuMP.is_valid(model::LinearFractionalModel, vref::VariableRef)
    return (model === vref.model &&
            vref.idx in keys(model.variables))
end

@forward LinearFractionalModel.model JuMP.num_variables

# Internal functions
#variable_info(vref::MyVariableRef) = vref.model.variables[vref.idx].info
# function update_variable_info(vref::MyVariableRef, info::JuMP.VariableInfo)
#     vref.model.variables[vref.idx] = JuMP.ScalarVariable(info)
# end

# JuMP.has_lower_bound(vref::MyVariableRef) = variable_info(vref).has_lb
# function JuMP.lower_bound(vref::MyVariableRef)::Float64
#     @assert !JuMP.is_fixed(vref)
#     return variable_info(vref).lower_bound
# end
# function JuMP.set_lower_bound(vref::MyVariableRef, lower)
#     info = variable_info(vref)
#     update_variable_info(vref,
#                          JuMP.VariableInfo(true, lower,
#                                            info.has_ub, info.upper_bound,
#                                            info.has_fix, info.fixed_value,
#                                            info.has_start, info.start,
#                                            info.binary, info.integer))
# end
# function JuMP.delete_lower_bound(vref::MyVariableRef)
#     info = variable_info(vref)
#     update_variable_info(vref,
#                          JuMP.VariableInfo(false, info.lower_bound,
#                                            info.has_ub, info.upper_bound,
#                                            info.has_fix, info.fixed_value,
#                                            info.has_start, info.start,
#                                            info.binary, info.integer))
# end
# JuMP.has_upper_bound(vref::MyVariableRef) = variable_info(vref).has_ub
# function JuMP.upper_bound(vref::MyVariableRef)::Float64
#     @assert !JuMP.is_fixed(vref)
#     return variable_info(vref).upper_bound
# end
# function JuMP.set_upper_bound(vref::MyVariableRef, upper)
#     info = variable_info(vref)
#     update_variable_info(vref,
#                          JuMP.VariableInfo(info.has_lb, info.lower_bound,
#                                            true, upper,
#                                            info.has_fix, info.fixed_value,
#                                            info.has_start, info.start,
#                                            info.binary, info.integer))
# end
# function JuMP.delete_upper_bound(vref::MyVariableRef)
#     info = variable_info(vref)
#     update_variable_info(vref,
#                          JuMP.VariableInfo(info.has_lb, info.lower_bound,
#                                            false, info.upper_bound,
#                                            info.has_fix, info.fixed_value,
#                                            info.has_start, info.start,
#                                            info.binary, info.integer))
# end
# JuMP.is_fixed(vref::MyVariableRef) = variable_info(vref).has_fix
# function JuMP.fix_value(vref::MyVariableRef)::Float64
#     return variable_info(vref).fixed_value
# end
# function JuMP.fix(vref::MyVariableRef, value; force::Bool = false)
#     info = variable_info(vref)
#     if !force && (info.has_lb || info.has_ub)
#         error("Unable to fix $(vref) to $(value) because it has existing bounds.")
#     end
#     update_variable_info(vref, JuMP.VariableInfo(
#         false, info.lower_bound, false, info.upper_bound, true, value,
#         info.has_start, info.start, info.binary, info.integer)
#     )
#     return
# end
# function JuMP.unfix(vref::MyVariableRef)
#     info = variable_info(vref)
#     update_variable_info(vref,
#                          JuMP.VariableInfo(info.has_lb, info.lower_bound,
#                                            info.has_ub, info.upper_bound,
#                                            false, info.fixed_value,
#                                            info.has_start, info.start,
#                                            info.binary, info.integer))
# end
# function JuMP.start_value(vref::MyVariableRef)::Union{Nothing, Float64}
#     return variable_info(vref).start
# end
# function JuMP.set_start_value(vref::MyVariableRef, start)
#     info = variable_info(vref)
#     update_variable_info(vref,
#                          JuMP.VariableInfo(info.has_lb, info.lower_bound,
#                                            info.has_ub, info.upper_bound,
#                                            info.has_fix, info.fixed_value,
#                                            true, start,
#                                            info.binary, info.integer))
# end
# JuMP.is_binary(vref::MyVariableRef) = variable_info(vref).binary
# function JuMP.set_binary(vref::MyVariableRef)
#     @assert !JuMP.is_integer(vref)
#     info = variable_info(vref)
#     update_variable_info(vref,
#                          JuMP.VariableInfo(info.has_lb, info.lower_bound,
#                                            info.has_ub, info.upper_bound,
#                                            info.has_fix, info.fixed_value,
#                                            info.has_start, info.start,
#                                            true, info.integer))
# end
# function JuMP.unset_binary(vref::MyVariableRef)
#     info = variable_info(vref)
#     update_variable_info(vref,
#                          JuMP.VariableInfo(info.has_lb, info.lower_bound,
#                                            info.has_ub, info.upper_bound,
#                                            info.has_fix, info.fixed_value,
#                                            info.has_start, info.start,
#                                            false, info.integer))
# end
# JuMP.is_integer(vref::MyVariableRef) = variable_info(vref).integer
# function JuMP.set_integer(vref::MyVariableRef)
#     @assert !JuMP.is_binary(vref)
#     info = variable_info(vref)
#     update_variable_info(vref,
#                          JuMP.VariableInfo(info.has_lb, info.lower_bound,
#                                            info.has_ub, info.upper_bound,
#                                            info.has_fix, info.fixed_value,
#                                            info.has_start, info.start,
#                                            info.binary, true))
# end
# function JuMP.unset_integer(vref::MyVariableRef)
#     info = variable_info(vref)
#     update_variable_info(vref,
#                          JuMP.VariableInfo(info.has_lb, info.lower_bound,
#                                            info.has_ub, info.upper_bound,
#                                            info.has_fix, info.fixed_value,
#                                            info.has_start, info.start,
#                                            info.binary, false))
# end

# Constraints
#const MyConstraintRef = JuMP.ConstraintRef{LinearFractionalModel, ConstraintIndex}

# TODO specialize on ScalarConstraint and VectorConstraint?

#
# function set_standard_form_coefficient(
#     constraint::ConstraintRef{Model, _MOICON{F, S}}, variable, value
#     ) where {S, T, F <: Union{MOI.ScalarAffineFunction{T}, MOI.ScalarQuadraticFunction{T}}}
#     MOI.modify(backend(owner_model(constraint)), index(constraint),
#                MOI.ScalarCoefficientChange(index(variable), convert(T, value)))
#     return
# end
# function set_standard_form_rhs(
#     constraint::ConstraintRef{Model, _MOICON{F, S}}, value) where {
#         T,
#         S <: Union{MOI.LessThan{T}, MOI.GreaterThan{T}, MOI.EqualTo{T}},
#         F <: Union{MOI.ScalarAffineFunction{T}, MOI.ScalarQuadraticFunction{T}}}
#     MOI.set(owner_model(constraint), MOI.ConstraintSet(), constraint,
#             S(convert(T, value)))
#     return
# end

# function transform_constraint(model::LinearFractionalModel, constraint::AbstractConstraint)
#     # Want to use set_coefficient, but this operates on a ConstraintRef (JuMP),
#     #, whereas add_constraint passes through an AbstractConstraint
#     # This is because we don't yet have the ConstraintRef since the model has not
#     # yet been added to JuMP.
#     #
#     # What if we instead
#     # 1.  Add the constraint to the JuMP Model
#     # 2.  Modify the resulting ConstraintRef
#
#     # Get the the scalar in standard form RHS
#
#
#
#     moi_function(c)
#
#
#
#     MOI.modify(backend(owner_model(constraint)), index(constraint),
#                MOI.ScalarCoefficientChange(index(variable), convert(T, value)))
# end

# TODO: remove after https://github.com/JuliaOpt/JuMP.jl/pull/1935 released
function set_standard_form_rhs(
    constraint::ConstraintRef{JuMP.Model, JuMP._MOICON{F, S}}, value) where {
        T,
        S <: Union{MOI.LessThan{T}, MOI.GreaterThan{T}, MOI.EqualTo{T}},
        F <: Union{MOI.ScalarAffineFunction{T}, MOI.ScalarQuadraticFunction{T}}}
    MOI.set(JuMP.owner_model(constraint), MOI.ConstraintSet(), constraint,
            S(convert(T, value)))
    return
end


# TODO: remove after https://github.com/JuliaOpt/JuMP.jl/pull/1935 released
function standard_form_rhs(
    constraint::ConstraintRef{JuMP.Model, JuMP._MOICON{F, S}}) where {
        T,
        S <: Union{MOI.LessThan{T}, MOI.GreaterThan{T}, MOI.EqualTo{T}},
        F <: Union{MOI.ScalarAffineFunction{T}, MOI.ScalarQuadraticFunction{T}}}
    con = JuMP.constraint_object(constraint)
    return MOIU.getconstant(con.set)
end

function transform_constraint(model::LinearFractionalModel, constraint_ref::ConstraintRef)
    α = standard_form_rhs(constraint_ref)
    set_standard_form_rhs(constraint_ref, 0.0)
    JuMP.set_coefficient(constraint_ref, model.t, -α)
end

JuMP.constraint_type(::LinearFractionalModel) = ConstraintRef
function JuMP.add_constraint(model::LinearFractionalModel, c::JuMP.AbstractConstraint,
                             name::String="")


    cref = JuMP.add_constraint(model.model, c, name)
    transform_constraint(model, cref)
    return cref
end

@forward LinearFractionalModel.model JuMP.delete, JuMP.is_valid
# function JuMP.delete(model::LinearFractionalModel, constraint_ref::MyConstraintRef)
#     @assert JuMP.is_valid(model.model, constraint_ref)
#     delete!(model.model.constraints, constraint_ref.index)
#     delete!(model.model.con_to_name, constraint_ref.index)
# end
# function JuMP.is_valid(model::LinearFractionalModel, constraint_ref::MyConstraintRef)
#     return (model.model === constraint_ref.model &&
#             constraint_ref.index in keys(model.model.constraints))
# end
# function JuMP.constraint_object(cref::MyConstraintRef)
#     return cref.model.constraints[cref.index]
# end

# Objective

function transform_objective(m::LinearFractionalModel)
    obj = objective_function(m.model, AffExpr)
    β = obj.constant
    obj.constant = 0.0
    JuMP.add_to_expression!(obj, β, m.t)
    JuMP.set_objective_function(m.model, obj)
end

function JuMP.set_objective(m::LinearFractionalModel, sense::MOI.OptimizationSense, f)
    JuMP.set_objective(m.model, sense, f)
    transform_objective(m)
end

macro numerator(model, args...)
    _error(str...) = _macro_error(:numerator, (model, args...), str...)

    # We don't overwrite `model` as it is used in `_error`
    esc_model = esc(model)
    if length(args) != 2
        # Either just an objective sense, or just an expression.
        _error("needs three arguments: model, objective sense (Max or Min) and expression.")
    end
    sense, x = args
    sense_expr = JuMP._moi_sense(_error, sense)
    newaff, parsecode = JuMP._parse_expr_toplevel(x, :q)
    code = quote
        q = Val{false}()
        $parsecode
        set_objective($esc_model, $sense_expr, $newaff)
    end
    return code
    #return JuMP._assert_valid_model(esc_internal_model, JuMP._macro_return(code, newaff))
end

macro denominator(model, expr)#args...)
    #_error(str...) = _macro_error(:objective, (model, args...), str...)

    # We don't overwrite `model` as it is used in `_error`
    esc_model = esc(model)
    # if length(args) != 1
    #     # Either just an objective sense, or just an expression.
    #     _error("needs two arguments: model and expression.")
    # end
    # sense, x = args
    # sense_expr = _moi_sense(_error, sense)
    newaff, parsecode = JuMP._parse_expr_toplevel(expr, :q)
    code = quote
        q = Val{false}()
        $parsecode
        set_denominator($esc_model, $newaff)
    end
    return code
    #return JuMP._assert_valid_model(esc_model.model, _macro_return(code, newaff))
end

function set_denominator(model::LinearFractionalModel, expr)
    if !isnothing(model.denominator_constraint)
        delete(model, model.denominator_constraint)
    end
    β = expr.constant
    expr.constant = 0.0
    add_to_expression!(expr, β, model.t)
    model.denominator_constraint = @constraint(model.model, expr == 1.0)
end

@forward LinearFractionalModel.model JuMP.objective_sense, JuMP.set_objective_sense
# JuMP.objective_sense(model::LinearFractionalModel) = JuMP.objective_sense(model.model.objectivesense)
# function JuMP.set_objective_sense(model::LinearFractionalModel, sense)
#     set_objective_sense(model.model, sense)
# end

@forward LinearFractionalModel.model JuMP.objective_function_type, JuMP.objective_function
#JuMP.objective_function_type(model::LinearFractionalModel) = JuMP.objective_function_type(model.model)
#JuMP.objective_function(model::LinearFractionalModel) = model.model.objective_function
# function JuMP.objective_function(model::LinearFractionalModel, FT::Type)
#     # InexactError should be thrown, this is needed in `objective.jl`
#     if !(model.model..objective_function isa FT)
#         throw(InexactError(:objective_function, FT,
#                            typeof(model.model.objective_function)))
#     end
#     return model.model.objective_function::FT
# end

# Names
# JuMP.name(vref::MyVariableRef) = vref.model.var_to_name[vref.idx]
# function JuMP.set_name(vref::MyVariableRef, name::String)
#     vref.model.var_to_name[vref.idx] = name
#     vref.model.name_to_var = nothing
# end

@forward LinearFractionalModel.model JuMP.variable_by_name

# function JuMP.variable_by_name(model::LinearFractionalModel, name::String)
#     if model.name_to_var === nothing
#         # Inspired from MOI/src/Utilities/model.jl
#         model.name_to_var = Dict{String, Int}()
#         for (var, var_name) in model.var_to_name
#             if haskey(model.name_to_var, var_name)
#                 # -1 is a special value that means this string does not map to
#                 # a unique variable name.
#                 model.name_to_var[var_name] = -1
#             else
#                 model.name_to_var[var_name] = var
#             end
#         end
#     end
#     index = get(model.name_to_var, name, nothing)
#     if index isa Nothing
#         return nothing
#     elseif index == -1
#         error("Multiple variables have the name $name.")
#     else
#         return MyVariableRef(model, index)
#     end
# end
# JuMP.name(cref::MyConstraintRef) = cref.model.con_to_name[cref.index]
# function JuMP.set_name(cref::MyConstraintRef, name::String)
#     cref.model.con_to_name[cref.index] = name
#     cref.model.name_to_con = nothing
# end

@forward LinearFractionalModel.model JuMP.constraint_by_name
# function JuMP.constraint_by_name(model::LinearFractionalModel, name::String)
#     if model.name_to_con === nothing
#         # Inspired from MOI/src/Utilities/model.jl
#         model.name_to_con = Dict{String, ConstraintIndex}()
#         for (con, con_name) in model.con_to_name
#             if haskey(model.name_to_con, con_name)
#                 # -1 is a special value that means this string does not map to
#                 # a unique constraint name.
#                 model.name_to_con[con_name] = ConstraintIndex(-1)
#             else
#                 model.name_to_con[con_name] = con
#             end
#         end
#     end
#     index = get(model.name_to_con, name, nothing)
#     if index isa Nothing
#         return nothing
#     elseif index.value == -1
#         error("Multiple constraints have the name $name.")
#     else
#         # We have no information on whether this is a vector constraint
#         # or a scalar constraint
#         return JuMP.ConstraintRef(model, index, JuMP.ScalarShape())
#     end
# end

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

@forward LinearFractionalModel.model optimize!

# Can't specialize `value` for LinearFractionalModel since it doesn't dispatch on the model-type
# and we re-use the variable_type.  Try at the MOI level
function value(v::VariableRef)::Float64
    return MOI.get(owner_model(v), MOI.VariablePrimal(), v)
end

# This doesn't work because `owner_model(v)` returns the inner model, not the LinearFractionalModel
function MOI.get(m::LinearFractionalModel, MOI.VariablePrimal(), v)
    MOI.get(m.model, MOI.VariablePrimal(), v)/MOI.get(m.model, MOI.VariablePrimal(), m.t)
end


end
