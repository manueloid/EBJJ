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

function orthogonality_check(n::Int64, η::ComplexF64, k::Float64)
    f(z::Float64) =
        time_dependent(n, η, k) * # this is the composition of the two time dependent parts of the product of the wave functions
        spatial_fourier(n, η, z) * # this is the left hand side of the integral
        ground_state(η, z) # this is the right hand side of the integral
    return quadgk(f, -1.0e3, 1.0e3, atol=1e-9)[1]
end

@testset "Testing ground state" begin
    zrange = range(-0.1, 0.1, length=1000)
    η = rand(ComplexF64)
    totest = [spatial_fourier(0, η, z) |> conj for z in zrange]
    reference = [ground_state(η, z) for z in zrange]
    @test isapprox(totest, reference, atol=1e-3)
end

@testset "orthogonality check for the general correction functions" begin
    # η = rand(ComplexF64)
    # k = rand(Float64)
    η = 1.0 + 0.5im
    k = 1.0
    for n in 0:5
        if n == 0
            @test isapprox(orthogonality_check(n, η, k), 1.0 + 0.0, atol=1e-3)
        else
            @test isapprox(orthogonality_check(n, η, k), 0.0 + 0.0, atol=1e-3)
        end
    end
end

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

function tdip_param(c::Control)
    ξ0, U = EBJJ.scaling_ξ0(c), c.U # constants
    b(t) = auxiliary(t, c)
    db(t) = ForwardDiff.derivative(b, t)
    itp = EBJJ.interpolation_integral(c)
    itpf = t::Float64 -> itp(t)
    η(t::Float64) = (ξ0^2 / b(t)^2 - 2im * db(t) / (U * b(t)))
    return η, itpf
end

function orthogonality_check(n::Int64, c::Control)
    η, itp = tdip_param(c)
    f(var) = time_dependent(n, conj(η(var[1])), itp(var[1])) * spatial_fourier(n, conj(η(var[1])), var[2]) * ground_state(η(var[1]), var[2])
    return hcubature(f, [0.0, -2.0e2], [c.T, 2.0e2], rtol=1e-3)[1]
end

@testset "orthogonality checks with system parameter" begin
    c = ControlFull()
    for n in 0:0
        if n == 0
            @test isapprox(orthogonality_check(n, c), c.T + 0.0, atol=1e-3)
        else
            @test isapprox(orthogonality_check(n, c), 0.0 + 0.0, atol=1e-3)
        end
    end
end
