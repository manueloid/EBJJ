using QuadGK
using Interpolations
using BenchmarkTools
using ForwardDiff
using SpecialPolynomials
using HCubature
using EBJJ

#=
# Orthogonality and normalisation checks

## 1 Orthogonality test for the general functions

in this document I will check what is the behaviour of the functions defined in [this file](./src/wavefunction.jl).

I defined the functions that will go in the integrand for the eSTA corrections, these functions can be used in such a way that it is possible to define the following composite function
But in this case, the dependency from the time is not explicit, and a complex variable is passed instead.
$$
    \chi^{*}_{n}(z,t) \chi_{0}(z,t)
$$

In the tests I want first to check the orthogonality for general input values and then I will pass the actual values from the system.

## 1 Orthogonality test

The way the functions are defined, I will be able to check either the normalisation of the ground state as well as the orthogonality of the ground state with respect to the excited states.
This is what I will do next.

### 1.1 Code implementatione

I will compose all the functions that I need and then integrate them numerically.
=#
#=
## 2 Orthogonality test for the actual system

Here I am going to pass the actual values of the system to the functions and check the orthogonality of the ground state with respect to the excited states.
By defining this function I will be able to make a blueprint for the actual integrals that I will have to evaluate in the future.

I will integrate also all over the time interval, the plan is to check that the normalisation of the ground state over the time interval is equal to the width of the time interval itself.

### 2.1 Code implementation

In this case I will define a function that will take the control parameter and the excitation number as the input and it will internally calculate all the relevant features.
The main features the function needs to calculate intenally are the following:
- The auxiliary function $b(t)$ and its derivative with respect to time
- The function which interpolates the integral function that goes into the time dependent phase of the wave function
- The function $ \eta (t) $ which is the argument of the Gaussian term of the STA wave function in momentum space

As I said, I can define everything inside the function call and then integrate numerically.
Actually, it would make sense to define a function that returns all the relevant features and then pass them to the integrand function.
I will call it `tdip_param` as in time dependent parameters.
=#

"""
    td_norm_test(c::Control, n::Int64)
Test the normalisation of the ground state for different times, given the control parameters and the number of points in the time interval.
"""
function td_norm_test(c::Control, n::Int64)
    trange = range(0.0, c.T, length=n)
    b(t) = auxiliary(t, c)
    db(t) = ForwardDiff.derivative(b, t)
    ξ0, U = EBJJ.scaling_ξ0(c), c.U # constants
    k(t::Float64) = EBJJ.interpolation_integral(c)(t)
    η(t::Float64) = (ξ0^2 / b(t)^2 - 2im * db(t) / (U * b(t)))
    f(z::Float64, t::Float64) = spatial_fourier(0, z, η(t)) * time_dependent(0, η(t), k(t)) * ground_state(z, η(t))
    return [quadgk(z -> f(z, t), -Inf, Inf, atol=1e-4)[1] for t in trange]
end
@testset "Testing normalisation for general `t`" begin
    c = ControlFull()
    n = 10
    result = td_norm_test(c, n)
    @test isapprox(result, ones(n), atol=1e-3)
end
"""
    td_norm_test(c::Control)
Test the normalisation of the ground state by integrating the two dimensional integral over the time variable as well as the position variable.
The plan is to check if the result is equal to the width of the time interval.
"""
function td_norm_test(c::Control)
    # trange = range(0.0, c.T, length=n)
    b(t) = auxiliary(t, c)
    db(t) = ForwardDiff.derivative(b, t)
    ξ0, U = EBJJ.scaling_ξ0(c), c.U # constants
    k(t::Float64) = EBJJ.interpolation_integral(c)(t)
    η(t::Float64) = (ξ0^2 / b(t)^2 - 2im * db(t) / (U * b(t)))
    f(v) = spatial_fourier(0, v[2], η(v[1])) * time_dependent(0, η(v[1]), k(v[1])) * ground_state(v[2], η(v[1]))
    result::ComplexF64, err::Float64 = hcubature(f, [0.0, -1.e3], [c.T, 1.e3], atol=1e-6)
    return result
end
@testset "Integrating all over the time interval" begin
    nparticles = 10:30:50
    tfs = 0.1pi:0.1pi:0.5pi
    for np in nparticles, tf in tfs
        c = ControlFull(np, tf)
        result = td_norm_test(c)
        @test isapprox(result, c.T, atol=1e-3)
    end
end

