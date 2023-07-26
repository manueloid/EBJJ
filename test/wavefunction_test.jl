#=
# Wave function testing

## 1 Testing of the phase integral function

In order to try to minimize the time it takes to calculate the corrections, I am going to interpolate the integral function and then define a new one.
Then I will test if the values of the functions are similar up to a certain tolerance.
=#

using QuadGK
using Interpolations
using BenchmarkTools

c = ControlFull(10, 0.2)
@btime EBJJ.imaginary_phase(2, 0.1, c)
