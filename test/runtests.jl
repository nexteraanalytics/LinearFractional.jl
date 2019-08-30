using Test
using LinearFractional
using LinearAlgebra
using Clp, Cbc
using JuMP

mytests = ["simpletest.jl", "arrayvars.jl", "mip.jl"]

for mytest in mytests
    include(mytest)
end
