# LinearFractional

<!-- [![Build Status](https://img.shields.io/travis/pluskid/Mocha.jl.svg?style=flat&branch=master)](https://travis-ci.org/pluskid/Mocha.jl)
[![Documentation Status](https://readthedocs.org/projects/mochajl/badge/?version=latest)](http://mochajl.readthedocs.org/)
[![Mocha](http://pkg.julialang.org/badges/Mocha_0.6.svg)](http://pkg.julialang.org/?pkg=Mocha&ver=0.6)
[![Coverage Status](https://img.shields.io/coveralls/pluskid/Mocha.jl.svg?style=flat)](https://coveralls.io/r/pluskid/Mocha.jl?branch=master) -->
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
<!--[![Build status](https://ci.appveyor.com/api/projects/status/342vcj5lj2jyegsp?svg=true)](https://ci.appveyor.com/project/pluskid/mocha-jl)-->

LinearFractional is an extension for [JuMP](https://github.com/JuliaOpt/JuMP.jl) to optimize linear programs with fractional objectives.  LinearFractional implements the [Charnes-Cooper transformation](https://en.wikipedia.org/wiki/Linear-fractional_programming) behind-the-scenes so that the user only needs to specify the problem as any ordinary JuMP problem, but specifying a numerator and denominator instead of a single objective function.


## Installation

To install the latest development version, run the following command:

```julia
Pkg.clone("PATH_TO_GIT_REPO")
```

Then you can run the built-in unit tests with

```julia
Pkg.test("LinearFractional")
```

to verify that everything is functioning properly on your machine.

## Basic Example

This toy example refers to the reference problem in http://www.ams.jhu.edu/~castello/625.414/Handouts/FractionalProg.pdf.

```julia
using LinearFractional
using Clp

lfp = LinearFractionalModel(solver=ClpSolver())
x1 = @variable(lfp, basename="x1", lowerbound=0)
x2 = @variable(lfp, basename="x2", lowerbound=0, upperbound=6)
@constraint(lfp, -x1 + x2 <= 4)
@constraint(lfp, 2x1 + x2 <= 14)
@constraint(lfp, x2 <= 6)
@numerator(lfp,  :Min, -2x1 + x2 + 2)
@denominator(lfp,  x1 + 3x2 + 4)
solve(lfp)
getobjectivevalue(lfp)
getvalue(x1)
getvalue(x2)
```

<!-- ## Documentation

The Mocha documentation is hosted at [readthedocs.org](http://mochajl.readthedocs.org/). -->
