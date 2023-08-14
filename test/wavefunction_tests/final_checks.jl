using QuadGK
using Interpolations
using BenchmarkTools
using ForwardDiff
using SpecialPolynomials
using HCubature
using EBJJ

#=
# Orthogonality and normalisation of the general form

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
        time_dependent(n, conj(η), k) * # this is the composition of the two time dependent parts of the product of the wave functions
        spatial_fourier(n, conj(η), z) * # this is the left hand side of the integral
        ground_state(η, z) # this is the right hand side of the integral
    return quadgk(f, -1.0e2, 1.0e2, rtol=1e-3)[1]
end

@testset "orthogonality check for the general correction functions" begin
    η = 1.0 + 0.0im
    k = 0.0
    for n in 0:5
        if n == 0
            @test isapprox(orthogonality_check(n, η, k), 1.0 + 0.0, atol=1e-3)
        else
            @test isapprox(orthogonality_check(n, η, k), 0.0 + 0.0, atol=1e-3)
        end
    end
end

#=
## Orthogonality test for the actual system

Here I am going to pass the actual values of the system to the functions and check the orthogonality of the ground state with respect to the excited states.
By defining this function I will be able to make a blueprint for the actual integrals that I will have to evaluate in the future.

I will integrate also all over the time interval, the plan is to check that the normalisation of the ground state over the time interval is equal to the width of the time interval itself.
=#
