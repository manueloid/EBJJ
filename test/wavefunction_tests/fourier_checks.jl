using QuadGK
using Interpolations
using BenchmarkTools
using ForwardDiff
using SpecialPolynomials
using HCubature
using EBJJ
using Test

#=
# Fourier checks
In this set of tests I am going to numerically check if the simplifications I made to turn the STA wave function from momentum to position representation are correct.

I have already checked in [the relevant repository](~/Repos/ExternalBJJ/fourier_transform.wlnb) that the Fourier transform of the product of a Gaussian function and of a Hermite polynomial has a specific close form.

Here I only want to check if the normalisation and the orthogonality hold.
This can be done by considering the fact that the Fourier transform of the whole STA wave function in momentum representation is the product of the Fourier transform of the spatial part and of the time dependent part.

Hence the STA wave function in position representation is just the product of the time dependent term and the Fourier transform of term dependent on the momentum varible $p$.

=#
#=
### 1 Fourier transform general complex number
I want to start from the general case, where the parameter $ \eta $ is a complex number and it is not tied to the specific of the system.
For the moment I do not want to be particularly efficient, I just want to check that the normalisation and the orthogonality hold.

I will thus define the general STA wave function and then test the normalisation for different $ n $.
I an not going to include the time dependent term, as it is not relevant for the normalisation.

=#

normalisation(η::ComplexF64, n::Int64) = (real(η) / pi)^(1 / 4) / sqrt(2^n * factorial(n)) * im^n * sqrt(conj(η) / η)^n / sqrt(η)
normalisation(n::Int64, η::ComplexF64) = normalisation(η, n)
gaussian(z::Float64, η::ComplexF64) = exp(-z^2 / (2 * η))
gaussian(η::ComplexF64, z::Float64) = gaussian(z, η)
hermite(n::Int64, z::Float64, η::ComplexF64) = SpecialPolynomials.basis(Hermite, n)(z * sqrt(real(η)) / abs(η))
hermite(n::Int64, η::ComplexF64, z::Float64) = hermite(n, z, η)
hermite(η::ComplexF64, n::Int64, z::Float64) = hermite(n, z, η)
spatial_fourier(n::Int64, z::Float64, η::ComplexF64) = normalisation(η, n) * gaussian(z, η) * hermite(n, z, η)
spatial_fourier(n::Int64, η::ComplexF64, z::Float64) = spatial_fourier(n, z, η)
spatial_fourier(η::ComplexF64, n::Int64, z::Float64) = spatial_fourier(n, z, η)

@testset "normalisation STA position general" begin
    η = rand(ComplexF64)
    for n in 0:5
        @test isapprox(quadgk(x -> conj(spatial_fourier(n, x, η)) * spatial_fourier(n, x, η), -Inf, Inf, atol=1e-7)[1], 1.0, atol=1e-3)
    end
end

#=
### 2 Fourier transform actual values 

Once the tests have been passed, I am going to use the same code but in this case I am going to pass the actual value of the parameter $ \eta $ at time $ t $.
In particular, we have that:
$$ \eta^2 = \left(- \frac{\xi_0^2}{b} + \frac{2i\dot{b}}{Ub} \right)$$

The plan is to test the normalisation.

#### 2.1 Normalisation of the full wave function

First I will check the orthogonality as it easier to do so, if that goes well I will continue with the normalisation, and implement multiple checks for it.

The only thing I need to do now is to define all the constant that appear in the STA wave function and then define some kind of $ \eta(t)$ and pass it to the `analytic` function.

All the relevant simplifications have been carried out in [the relevant notebook](~/Repos/ExternalBJJ/wrapup.wlnb), so there is no need to write them again here.
=#

"""
    `spatial_fourier(n, t, z, c::Control)`
Return the value of the STA wave function in position representation at the point `x` and at time `t`, given the control parameters `c`.
This is the function with absolutely no simplifications and it is the one I am going to test against.
There is no time dependent part in this function, as it is not relevant for the normalisation.
"""
function wave_function(n::Int64, t, x, c::Control)
    J0, N, U = c.J0, c.N, c.U
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    α(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) - im * b(t) / (U * db(t)) # Parameter of the Gaussian term
    return spatial_fourier(n, x, α(t))
end
@testset "normalisation with system parameters" begin
    c = ControlFull(10, 0.03)
    norm(n, t) = quadgk(x -> conj(wave_function(n, t, x, c)) * wave_function(n, t, x, c), -Inf, Inf, atol=1e-7)[1] |> real
    for t in range(1e-3, c.T - 1e-3, length=100), n in 0:5
        @test isapprox(norm(n, t), 1.0, atol=1e-3)
    end
end

#=
### 3 Splitting all the stuff
Here I need to define the functions that will go in the calculation of the corrections.
I basically have something of the form $ \langle \chi_n| \hat{O} | \chi_0  \rangle $, where $ \hat{O} $ is some kind of operator.

The righthand side of the integral is nothing more that the ground state, while for the lefthand side I need to take the complex conjugate of the general STA wave function.

Code-wise, I think I am going to define a function to evaluate the lefthand side where I will pass the complex conjugate of the variable instead of taking the complex conjugate of the whole function.

First let me check if the two approaches are equivalent.
=#

function wave_function_conj(n::Int64, t, x, c::Control)
    J0, N, U = c.J0, c.N, c.U
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    α(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) + im * b(t) / (U * db(t)) # Parameter of the Gaussian term
    return spatial_fourier(n, x, α(t))
end
@testset "conjugation" begin
    c = ControlFull(10, 0.03)
    for _ in 1:10000
        x = rand() * 10
        t = c.T * rand()
        n = rand(0:2:10)
        @test isapprox(wave_function_conj(n, t, x, c), conj(wave_function(n, t, x, c)), atol=1e-7)
    end
end

#=
### 4 Lefthand side

Here I am going to set up the code for the left hand side of the integrand.
I will write down the calculations in a separate notebook, so here I will just write down the code.

First I will define everything in terms of a general complex variable and then I will pass the actual value of the parameter $ \eta $.
=#

lhs_normalisation(η::ComplexF64, n::Int64) = (real(η) / pi)^(1 / 2) / sqrt(2^n * factorial(n)) * (-im)^n * sqrt(conj(η) / η)^n / abs(η)
lhs_normalisation(n::Int64, η::ComplexF64) = lhs_normalisation(η, n)
lhs_spatial_fourier(n::Int64, z::Float64, η::ComplexF64) = lhs_normalisation(η, n) * gaussian(z, η) * hermite(n, z, η)
# When passing the actual value of the parameter η, I need to take the complex conjugate of the parameter
function test_ground_state(η::ComplexF64,n::Int64)
    lhs(x) = lhs_spatial_fourier(n, x, conj(η))
    rhs(x) = gaussian(x, η)
    return quadgk(x -> lhs(x) * rhs(x), -Inf, Inf,  atol=1e-7)[1]
end

@testset "Testing lhs general complex" begin
    for _ in 1:10000
        η = rand(ComplexF64)
        n = rand(0:2:10)
        if n == 0
            @test isapprox(test_ground_state(η, n), 1.0, atol=1e-3)
        else
            @test isapprox(test_ground_state(η, n), 0.0, atol=1e-3)
        end
    end
end

function test_ground_state(n::Int64, t, c::Control)
    J0, N, U = c.J0, c.N, c.U
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    α(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) - im * b(t) / (U * db(t)) # Parameter of the Gaussian term
    return test_ground_state(α(t), n)
end

@testset "Testing lhs actual values" begin
    c = ControlFull(10, 0.03)
    for t in range(1e-3, c.T - 1e-3, length=100), n in 0:5
        if n == 0
            @test isapprox(test_ground_state(n, t, c), 1.0, atol=1e-3)
        else
            @test isapprox(test_ground_state(n, t, c), 0.0, atol=1e-3)
        end
    end
end

