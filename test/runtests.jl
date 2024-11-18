using Test
using LinearFractional
using LinearAlgebra
using Clp, Cbc
using ECOS
using JuMP

mytests = ["simpletest.jl", "arrayvars.jl", "mip.jl", "socp.jl"]

for mytest in mytests
    include(mytest)
end
