using QuadGK
using Interpolations
using BenchmarkTools
using ForwardDiff
using SpecialPolynomials
using HCubature
using EBJJ

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
    norm(n::Int64, η::ComplexF64) = quadgk(z -> conj(spatial_fourier(n, z, η)) * spatial_fourier(n, z, η), -Inf, Inf, atol=1e-7)[1]
    for n in 0:5
        @test isapprox(norm(n, η), 1.0, atol=1e-7)
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
    `wave_function(n, t, x, c::Control)`
Return the value of the STA wave function in position representation at the point `x` and at time `t`, given the control parameters `c`.
This is the function with absolutely no simplifications and it is the one I am going to test against.
"""
function wave_function(n::Int64, t, x, c::Control)
    ξ0, U = EBJJ.scaling_ξ0(c), c.U # constants
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    normalisation = sqrt(EBJJ.scaling_ξ0(c)) / (π)^0.25 # Normalisation factor
    itp = EBJJ.interpolation_integral(c) # Interpolating function for the integral of the phase
    imaginary_phase(n, t) = exp(-im * (n + 0.5) * ξ0^2 * U / 2 * itp(t)) # Imaginary phase
    η(t) = (ξ0^2 / b(t)^2 - 2im * db(t) / (U * b(t))) # Argument of the Gaussian term
    spatial(n, t, x) = analytic(n, x, η(t)) # Spatial term
    excitation(n, t) = sqrt(2^n * factorial(n) * b(t)) # Excitation term
    return normalisation / excitation(n, t) * imaginary_phase(n, t) * spatial(n, t, x)
end

#=
#=
#### 2.1.1 Code implementation

I think it would make sense to integrate over the two dimensions $ x $ and $ t $, so that we can check the normalisation for different values of $ t $.

If the normalisation is correct, we should have that:

$$ \int_{0}^{t_f} \int_{-\infty}^{\infty} \psi_n^*(x, t) \psi_n(x, t) dx dt = t_f $$

where $ t_f $ is the final time of the STA protocol.
I am goinot to use low tolerance, as the requirement for the normalisation is not too strict
=#

@testset "normalisation STA position full" begin
    c = ControlFull()
    for n in 0:5
        # defining a new function that takes only two arguments
        f(x) = conj(wave_function(n, x[1], x[2], c)) * wave_function(n, x[1], x[2], c)
        # integrating over the two dimensions
        normalisation = hcubature(f, [0.0, -1e3], [c.T, 1e3], rtol=1e-3)[1]
        @test isapprox(normalisation |> real, c.T, atol=1e-3)
    end
end
=#


