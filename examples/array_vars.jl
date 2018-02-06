# Validated from http://www.ams.jhu.edu/~castello/625.414/Handouts/FractionalProg.pdf
#using LinearFractional
include(joinpath(Pkg.dir("LinearFractional")"/src/LinearFractional.jl"))
using Clp
using JuMP

lfp = LinearFractional.LinearFractionalModel(solver=ClpSolver())
x = @variable(lfp, [i=1:2], basename="x", lowerbound=0)
@constraint(lfp, x[2] <= 6)
@constraint(lfp, -x[1] + x[2] <= 4)
@constraint(lfp, 2x[1] + x[2] <= 14)
LinearFractional.@numerator(lfp,  :Min, -2x[1] + x[2] + 2)
LinearFractional.@denominator(lfp,  x[1] + 3x[2] + 4)
solve(lfp)
getobjectivevalue(lfp)
getvalue(x)
