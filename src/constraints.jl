struct LinearFractionalConstraint <: JuMP.AbstractConstraint
    lctrans::Vector{LinearConstraint}  # 1 or 2 (lb, ub)
end
const LinearFractionalConstraintRef = ConstraintRef{LinearFractionalModel, LinearFractionalConstraint}


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
    LinearFractionalConstraint(cons)
end


function addconstraint(model::LinearFractionalModel, constraint::LinearFractionalConstraint)
    for con in constraint.lctrans
        addconstraint(model.transformedmodel, con)
    end
    LinearFractionalConstraintRef(model, length(constraint.lctrans))
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

# _vectorize_like(x::Number, y::AbstractArray{LinearFractionalAffExpr}) = (ret = similar(y, typeof(x)); fill!(ret, x))
# function _vectorize_like(x::AbstractArray{R}, y::AbstractArray{LinearFractionalAffExpr}) where R<:Number
#     for i in 1:max(ndims(x),ndims(y))
#         _size(x,i) == _size(y,i) || error("Unequal sizes for ranged constraint")
#     end
#     x
# end
#
# function constructconstraint!(x::AbstractArray{LinearFractionalAffExpr}, lb, ub)
#     LB = _vectorize_like(lb,x)
#     UB = _vectorize_like(ub,x)
#     ret = similar(x, LinearConstraint)
#     map!(ret, eachindex(ret)) do i
#         constructconstraint!(x[i], LB[i], UB[i])
#     end
# end
