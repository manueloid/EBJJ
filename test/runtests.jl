using EBJJ
using Test
#= 
# Tests for the `EBJJ.jl` package


In the following I will test the functions in the `EBJJ.jl` package.
I will focus on the polynomial first and then I will move onto the wave function.
Hopefully, in the end, I will be able to test the corrections to understand if they are correct.

In the remainder, I will only focus on the ControlFull type for the moment.
Just remember that the default parameters are:

- $ N = 30 $
- $ T = 0.05\pi $
- $ J_0 = 0.1 $
- $ J_f = 0.05 $
- $ U = 0.02 $


=#
# include("control.jl")
include("./wavefunction_tests/momentum_checks.jl")
# include("./wavefunction_tests/fourier_checks.jl")
