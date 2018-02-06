using Base.Test
using LinearFractional
using Clp
using JuMP

mytests = ["simpletest.jl", "arrayvars.jl"]

for mytest in mytests
    include(mytest)
end
