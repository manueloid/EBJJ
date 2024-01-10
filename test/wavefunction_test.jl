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

using SpecialPolynomials, QuadGK, EBJJ, Test, HCubature
he(n::Int64, x::Float64) = SpecialPolynomials.basis(Hermite, n)(x)
"""
    momentum(n::Int64, p::Float64, α::ComplexF64)
Wave function in momentum representation, for a general complex number `α`.
This function does not have the imaginaray phase term.
"""
function momentum(n::Int64, z::Float64, α::ComplexF64, h::Float64)
    r = sqrt(real(α))
    f(p, z) = (2pi * h)^(-1 / 2) * exp(im * p * z / h) * # Fourier transform term
              (r^2 / pi)^(1 / 4) * # Normalisation term
              (2^n * factorial(n))^(-1 / 2) * # Hermite polynomial normalisation term
              exp(-p^2 * α / 2) * # Gaussian term
              he(n, p * r) # Hermite polynomial
    return quadgk(p -> f(p, z), -Inf, Inf, atol=1e-7)[1]
end
"""
    position(n::Int64, z::Float64, α::ComplexF64, h::Float64)
Fourier transform of the wave function in momentum representation, for a general complex number `α`.
In this case I needed to introduce the parameter `h` coming from the Fourier transform.
This function does not have the imaginaray phase term.
"""
function position(n::Int64, z::Float64, α::ComplexF64, h::Float64)
    r = sqrt(real(α))
    return (r^2 / pi)^(1 / 4) * # Normalisation term
           (2^n * factorial(n))^(-1 / 2) * # Hermite polynomial normalisation term
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
    φ(t::Float64, c::Control)
Function that returns the value of the phase factor for the given control parameters.
"""
function φ(t::Float64, c::Control)
    h, Λ = 2.0 / c.N, EBJJ.Λ(c)
    imag_phase_integrand(t::Float64) = 2Λ * h * real(f2(t, c)) # Integrand of the phase factor
    return quadgk(τ -> imag_phase_integrand(τ), 0.0, t, atol=1e-7)[1]
end
"""
    momentum(n::Int64, z::Float64, t::Float64, c::Control)
Function that returns the Fourier transform of the wave function in momentum representation for a certain time `t`, given the parameters of the system in the `Control` type.
"""
function momentum(n::Int64, z::Float64, t::Float64, c::Control)
    h = 2.0 / c.N
    α = f2(t, c)
    return momentum(n, z, α, h) * exp(-im * (n + 1 / 2) * φ(t, c)) # Phase factor
end
"""
    position(n::Int64, z::Float64, t::Float64, c::Control)
Function that returns the Fourier transform of the wave function in position representation for a certain time `t`, given the parameters of the system in the `Control` type.
"""
function position(n::Int64, z::Float64, t::Float64, c::Control)
    h = 2.0 / c.N
    α = f2(t, c)
    return position(n, z, α, h) * exp(-im * (n + 1 / 2) * φ(t, c)) # Phase factor
end

@testset "Testing the Fourier transform for f²" begin
    for _ in 1:5000
        n, z, t = rand(0:2:8), rand(), rand()
        N = rand(10:10:50)
        h = 2.0 / N
        c = ControlFull(N, rand(), rand(), t, 2, 2:2)
        @test isapprox(momentum(n, z, t, c), position(n, z / h, t, c), atol=1e-4)
    end
end

#=
### 2. Corrections 

In this section I will check if the simplifications I made give the same result as the full expression.
I will first test the corrections for a general complex number and then I will pass the actual value of the parameter $ f^2$..

The wave function that I am going to use - the $ \chi_n(z,t) $ - is the `position` function defined above.

The plan is to see if the integral over the whole 2D interval of $\chi_n(z,t) \chi_0(z,t)$ gives the same result as the integral over the whole 2D interval of the simplified version that I got elsewhere.

After that I will try to see if the value of the corrections are the same.

#### 2.1. Corrections for a general complex number

I have the `position` function defined above, that is the full wave function, so I only need to take the product of that function with the complex conjugate of another one for different `n`.
I will get the simplified version of the product of two wave functions from another file and then compare the two integrals.     
For this case I will not pass the imaginary phase integral to the `position` function, as I will need it later.
=#

"""
    simplified(n::Int64, α2::ComplexF64, α2c::ComplexF64, h::Float64)
Return the time dependent part of the product ⟨χₙ|χ₀⟩, for a general complex number `α` its complex conjugate `αc` and the energy level `n`.
This is the simplified version of the product of two wave functions.
The function also takes the parameter `h`, for later use.
The function takes two complex numbers as input to reduce the computational cost.
"""
function simplified(n::Int64, z::Float64, α::ComplexF64, h::Float64)
    r = sqrt(real(α))
    αc = conj(α)
    return (r^2 / pi)^(1 / 2) * # Normalisation term
           (2^n * factorial(n))^(-1 / 2) * # Hermite polynomial normalisation term
           ((-im)^n / h) / abs(α) * (α / αc)^(n / 2) * # Fourier transform normalisation term
           exp(-z^2 / (2 * αc)) * he(n, z * r / abs(α)) * exp(-z^2 / (2 * α))
end
"""
    non_simplified(n::Int64, z::Float64, α::ComplexF64,  h::Float64)
Return the time dependent part of the product ⟨χₙ|χ₀⟩, for a general complex number `α` and the energy level `n`, without any simplification.
"""
function non_simplified(n::Int64, z::Float64, α::ComplexF64, h::Float64)
    αc = conj(α)
    return position(n, z, αc, h) * position(0, z, α, h)
end

# Simple test to see if the two functions give the same result, without integration
@testset "Testing the simplified version" begin
    for _ in 1:500
        n, z, α, h = rand(0:2:8), rand(), rand() + im * rand(), rand()
        @test isapprox(simplified(n, z / h, α, h), non_simplified(n, z / h, α, h), atol=1e-4)
    end
end

#=
#### 2.2. Corrections for $ f^2$

Now I will do the same thing as before, but I will pass the actual value of the parameter $ f^2$ for a certain time `t`.
It is nothing different from what I did in section 1.2.
=#

"""
    simplified(n::Int64, z::Float64, t::Float64, c::Control)
Return the time dependent part of the product ⟨χₙ|χ₀⟩, for a certain time `t`, given the parameters of the system in the `Control` type.
This is the simplified version of the product of two wave functions.
"""
function simplified(n::Int64, z::Float64, t::Float64, c::Control)
    h = 2.0 / c.N
    α = f2(t, c)
    return simplified(n, z, α, h) * exp(-im * (n + 1 / 2) * φ(t, c)) # Phase factor
end
"""
    non_simplified(n::Int64, z::Float64, t::Float64, c::Control)
Return the time dependent part of the product ⟨χₙ|χ₀⟩, for a certain time `t`, given the parameters of the system in the `Control` type, without any simplification.
"""
function non_simplified(n::Int64, z::Float64, t::Float64, c::Control)
    h = 2.0 / c.N
    α = f2(t, c)
    return non_simplified(n, z, α, h) * exp(-im * (n + 1 / 2) * φ(t, c)) # Phase factor
end

@testset "Testing the simplified version for f²" begin
    for _ in 1:5000
        n, z, t = rand(0:2:8), rand(), rand()
        N = rand(10:10:50)
        h = 2.0 / N
        c = ControlFull(N, rand(), rand(), t, 2, 2:2)
        @test isapprox(simplified(n, z, t, c), non_simplified(n, z, t, c), atol=1e-4)
    end
end

#=
#### 2.3. Testing the values of the integrals

Now I will test if the integrals of the two functions give the same result.
I will do that directly for the value of $ f^2 $.

I will both perform a 2d integration and a double 1d integration.
=#

@testset "Testing the integrals for f²" begin
    for _ in 1:5
        n, t = rand(0:2:8), rand()
        N = rand(10:10:50)
        h = 2.0 / N
        c = ControlFull(N, rand(), rand(), t, 2, 2:2)
        simplified_integrand(t) = quadgk(z -> simplified(n, z, t, c), -1., 1., atol=1e-7)[1] # 2d integration
        non_simplified_integrand(t) = quadgk(z -> non_simplified(n, z, t, c), -1., 1., atol=1e-7)[1] # 2d integration
        @test isapprox(
            quadgk(simplified_integrand, 0.0, t, atol=1e-7)[1], # double 1d integration
            quadgk(non_simplified_integrand, 0.0, t, atol=1e-7)[1], # double 1d integration
            atol=1e-4
        )
        @test isapprox(
            hcubature(v -> simplified(n, v[1], v[2], c), [-10., 0.0], [1, t], atol=1e-7)[1], # double 1d integration
            hcubature(v -> non_simplified(n, v[1], v[2], c), [-10., 0.0], [1, t], atol=1e-7)[1], # double 1d integration
            atol=1e-4
        )
    end
end

#=
### 3. Conclusions

The tests are all passed, so I can start from here and use the simplified version of the wave function to calculate the corrections 

=#

