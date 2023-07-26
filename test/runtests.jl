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
@testset begin
    # Need to check if the gaussian_arg function and its conjugate work correctly
    c = ControlFull()
    α(t) = EBJJ.gaussian_arg(t, c)
    αc(t) = conj(α(t))
    @test real(α(0.1)) == real(αc(0.1))
    @test imag(α(0.1)) == -imag(αc(0.1))
end
