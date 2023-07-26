#=
# Main file for the EBJJ package.

In this file I will define the module EBJJ, by including all the files I deem necessary.
I then import the external packages needed for the package to work.
Finally, I will export the functions that I want to be available to the user.

Here is a short breakdown of the files that I will import:

- `base.jl` is the file where I define the parent type `Control` and the derived types `ControlFull` and `ControlInt`, plus some constructor to speed up the workflow

Then, for a list of the packages I need, here is the list and where I use them:
=#
module EBJJ

# Packages used in development
using Revise # This is just used during development
using Cthulhu # This is used to speed up the code
using BenchmarkTools # This is used to benchmark the code

# Packages used in the code
using ForwardDiff
using QuadGK
using Interpolations
using SpecialPolynomials

include("base.jl")
export Control
export ControlFull
export ControlInt

include("polynomials.jl")
export auxiliary
export control_function

include("wavefunction.jl")

end
