using Base.Test
using LinearFractional
using Clp
using JuMP

mytests = ["simpletest.jl"]

for mytest in mytests
    include(mytest)
end
