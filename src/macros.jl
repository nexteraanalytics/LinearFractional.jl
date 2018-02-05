# Returns a new variable belonging to the model `m`. Additional positional arguments can be used to dispatch the call to a different method.
# The return type should only depends on the positional arguments for `variabletype` to make sense.

variabletype(m::LinearFractionalModel) = LinearFractionalVariable


function constructvariable!(m::LinearFractionalModel, _error::Function, lowerbound::Number, upperbound::Number, category::Symbol, objective::Number, inconstraints::Vector, coefficients::Vector{Float64}, basename::AbstractString, start::Number; extra_kwargs...)
    for (kwarg, _) in extra_kwargs
        _error("Unrecognized keyword argument $kwarg")
    end
    LinearFractionalVariable(m, lowerbound, upperbound, category == :Default ? :Cont : category, objective, inconstraints, coefficients, basename, start)
end

function constructvariable!(m::LinearFractionalModel, _error::Function, lowerbound::Number, upperbound::Number, category::Symbol, basename::AbstractString, start::Number; extra_kwargs...)
    for (kwarg, _) in extra_kwargs
        _error("Unrecognized keyword argument $kwarg")
    end
    LinearFractionalVariable(m, lowerbound, upperbound, category == :Default ? :Cont : category, basename, start)
end


function constructconstraint!(aff::LinearFractionalAffExpr, sense::Symbol)
    #offset = aff.constant * aff.t
    #aff.constant = 0.0

    if sense == :(<=) || sense == :≤
        return LinearConstraint(aff, -Inf, 0)
    elseif sense == :(>=) || sense == :≥
        return LinearConstraint(aff, 0, Inf)
    elseif sense == :(==)
        return LinearConstraint(aff, 0, 0)
    else
        error("Cannot handle ranged constraint")
    end
end


function constructconstraint!(aff::LinearFractionalAffExpr, lb, ub)
#    LinearConstraint(aff, lb-offset, ub-offset)
    return LinearConstraint(aff, lb, ub)
    # LinearConstraint(aff - lb * aff.t, 0, Inf)
    # LinearConstraint(aff - ub * aff.t, -Inf, 0)
end


macro objective(m, args...)
    m = esc(m)
    if length(args) != 2
        # Either just an objective sene, or just an expression.
        error("in @objective: needs three arguments: model, objective sense (Max or Min) and expression.")
    end
    sense, x = args
    if sense == :Min || sense == :Max
        sense = Expr(:quote,sense)
    end
    newaff, parsecode = JuMP.parseExprToplevel(x, :q)
    code = quote
        q = Val{false}()
        $parsecode
        setobjective($m, $(esc(sense)), $newaff)
    end
    #return newaff, parsecode
    return assert_validmodel(m, code)
end

macro lfobjective(m, args...)
    m = esc(m)
    if length(args) != 3
        # Either just an objective sense, or just an expression.
        error("in @objective: needs three arguments: model, objective sense (Max or Min) and expression.")
    end
    sense, numer, denom = args
    if sense == :Min || sense == :Max
        sense = Expr(:quote, sense)
    end
    numeraff, numerparsecode = JuMP.parseExprToplevel(numer, :q)
    denomaff, denomparsecode = JuMP.parseExprToplevel(denom, :q)
    code = quote
        q = Val{false}()
        $numerparsecode
        #setobjective($m, $(esc(sense)), $numeraff, $denomaff)
        setobjective($(m).transformedmodel, $(esc(sense)), $(numeraff).afftrans)
    end

    #return code
    #return JuMP.assert_validmodel(m, code)
end

macro numerator(m, args...)
    m = esc(m)
    if length(args) != 2
        # Either just an objective sense, or just an expression.
        error("in @objective: needs three arguments: model, objective sense (Max or Min) and expression.")
    end
    sense, numer = args
    if sense == :Min || sense == :Max
        sense = Expr(:quote, sense)
    end
    numeraff, numerparsecode = JuMP.parseExprToplevel(numer, :q)
    code = quote
        q = Val{false}()
        $numerparsecode
        #setobjective($m, $(esc(sense)), $numeraff, $denomaff)
        setobjective($(m).transformedmodel, $(esc(sense)), $(numeraff).afftrans)
    end
    #return code
    return JuMP.assert_validmodel(m, code)
end

macro denominator(m, denom)
    m = esc(m)
    denomaff, denomparsecode = JuMP.parseExprToplevel(denom, :q)
    code = quote
        q = Val{false}()
        $denomparsecode
        #setobjective($m, $(esc(sense)), $numeraff, $denomaff)
        setdenominator!($(m), $(denomaff))
    end
    return code
    #return JuMP.assert_validmodel(m, code)
end


# function JuMP.constructvariable!(m::LinearFractionalModel, _error::Function, haslb::Bool, lowerbound::Number, hasub::Bool, upperbound::Number,
#                             hasfix::Bool, fixedvalue::Number, binary::Bool, integer::Bool, name::AbstractString,
#                             hasstart::Bool, start::Number; extra_kwargs...)
#     for (kwarg, _) in extra_kwargs
#          _error("Unrecognized keyword argument $kwarg")
#      end
#      transvar = Variable(m.transformedmodel)
#      v = LinearFractionalVariable(transvar, m)
#      if haslb
#          setlowerbound(v, lowerbound)
#      end
#      if hasub
#          setupperbound(v, upperbound)
#      end
#      # if hasfix
#      #     fix(v, fixedvalue)
#      # end
#      # if binary
#      #     setbinary(v)
#      # end
#      # if integer
#      #     setinteger(v)
#      # end
#      # TODO: MOIU.Instance does not support VariablePrimalStart
#      #if hasstart
#      #    setstartvalue(v, start)
#      #end
#      if name != EMPTYSTRING
#          setname(v, name)
#      end
#      return v
# end
