#=
# Testing the Wave function

In the following code, I will try to test if the wave function and all the simplifications are working properly.
I am not going to repeat here all the calculations, the notes are there for a reason.

I will first test the wave function for a general complex number and then I will pass the actual value of the parameter $ f^2$.

The first thing I need to test is the Fourier transform of the wave function.
=#
#=
### 1. Fourier transform of the wave function

#### 1.1. Fourier transform of the wave function for a general complex number

I am not going to use any simplifications here, I will just define the wave function in momentum representation and that will be the function I will pass to the `quadgk` function.
Furthermore, I am not going to split the wave function in smaller chunks, I will just define the full function, trying to split it on multiple lines to make it more readable.
The only thing I am going to define is the Hermite polynomial, as I will need it later.

I will define the function `momentum` which is the numerical Fourier transform of the wave function in momentum representation.
Then I will define the function `position` which is the analytic solution of the Fourier transform of the wave function in momentum representation.

I will finally test the two functions for random values of the parameters.
=#

using SpecialPolynomials, QuadGK, EBJJ, Test
he(n::Int64, x::Float64) = SpecialPolynomials.basis(Hermite, n)(x)
"""
    momentum(n::Int64, p::Float64, α::ComplexF64)
Wave function in momentum representation, for a general complex number `α`.
"""
function momentum(n::Int64, z::Float64, α::ComplexF64, h::Float64)
    r = sqrt(real(α))
    f(p, z) = (2pi * h)^(-1 / 2) * exp(im * p * z / h) * # Fourier transform term
              (r^2 / pi)^(1 / 4) * # Normalisation term
              (2^n * factorial(n))^(-1 / 2) * # Hermite polynomial normalisation term
              exp(-im * (n + 1 / 2) * atan(1 / r)) * # Phase term
              exp(-p^2 * α / 2) * # Gaussian term
              he(n, p * r) # Hermite polynomial
    return quadgk(p -> f(p, z), -Inf, Inf, atol=1e-7)[1]
end
"""
    position(n::Int64, z::Float64, α::ComplexF64, h::Float64)
Fourier transform of the wave function in momentum representation, for a general complex number `α`.
In this case I needed to introduce the parameter `h` coming from the Fourier transform.
"""
function position(n::Int64, z::Float64, α::ComplexF64, h::Float64)
    r = sqrt(real(α))
    return (r^2 / pi)^(1 / 4) * # Normalisation term
           (2^n * factorial(n))^(-1 / 2) * # Hermite polynomial normalisation term
           exp(-im * (n + 1 / 2) * atan(1 / r)) * # Phase term
           im^n / sqrt(h) * 1 / sqrt(α) * (conj(α) / α)^(n / 2) * # Fourier transform term
           exp(-z^2 / (2 * α)) * # Gaussian term
           he(n, z * r / abs(α)) # Hermite polynomial
end

@testset "Testing the Fourier transform for random values" begin
    for _ in 1:5000
        n, z, α, h = rand(0:2:8), rand(), rand() + im * rand(), rand()
        @test isapprox(momentum(n, z, α, h), position(n, z / h, α, h), atol=1e-4)
    end
end

#=
#### 1.2. Fourier transform for $ f^2$

Now I will do the same thing as before, but I will pass the actual value of the parameter $ f^2$ for a certain time `t`.

I need to define a function that returns the value of $ f^2$ for a certain time `t`, given the parameters of the simulation, stored in the `Control` type.
=#

"""
    f2(t::Float64, c::Control)
Return the value of the argument of the Gaussian term of the wave function for a certain time `t`, given the parameters of the system in the `Control` type.
"""
function f2(t::Float64, c::Control)
    h, Λ = 2.0 / c.N, EBJJ.Λ(c)
    # Definition of the auxiliary functions
    b(t) = auxiliary(t, c)
    db(t) = EBJJ.auxiliary_1d(t, c)
   return sqrt(1 / (2Λ)) * 1 / (h * b(t)^2) - im * db(t) / (2Λ * h * b(t))
end
"""
    momentum(n::Int64, z::Float64, t::Float64, c::Control)
Function that returns the Fourier transform of the wave function in momentum representation for a certain time `t`, given the parameters of the system in the `Control` type.
"""
function momentum(n::Int64, z::Float64, t::Float64, c::Control)
    h = 2.0 / c.N
    α = f2(t, c)
    return momentum(n, z, α, h)
end
"""
    position(n::Int64, z::Float64, t::Float64, c::Control)
Function that returns the Fourier transform of the wave function in position representation for a certain time `t`, given the parameters of the system in the `Control` type.
"""
function position(n::Int64, z::Float64, t::Float64, c::Control)
    h = 2.0 / c.N
    α = f2(t, c)
    return position(n, z / h, α, h)
end

@testset "Testing the Fourier transform for f²" begin
    for _ in 1:5000
        n, z, t = rand(0:2:8), rand(), rand() 
        c = ControlFull( rand(10:10:50), rand(), rand(), t, 2, 2:2)
        @test isapprox(momentum(n, z, t, c), position(n, z, t, c), atol=1e-4)
    end
end
