using QuadGK
using Interpolations
using BenchmarkTools
using ForwardDiff
using SpecialPolynomials
using EBJJ

#=
# Fourier checks
In this set of tests I am going to numerically check if the subsitutions I made to turn the STA wave function from momentum to position representation are correct.

My plan is to find make as much simplifications as possible, but in order to do that, I need to know if the substitutions I made are correct.

Plans:
I have already checked that the STA wf in monentum representation is correct, so I will use that as a reference.
I will hence perform a numerical Fourier transform on the STA wf in momentum representation and then I will compare it with the analytical solution that I got from the internet.

First, I will define the analytical solution of the Fourier transform of the product of a Gaussian function of the form $ e^{-p^2/2\alpha}$ and a hermite polynomial of the form $ \mathcal_n(\beta p) $

It is given by:
$$ i^n \sqrt{\alpha} \exp\left\{\frac { -\alpha x^2}{2} \right\} \gamma^n \mathcal{H}_n\left( \frac{\alpha\beta}{\gamma}x\right)$$
where
$ \gamma = \sqrt{2 \alpha \beta^2 - 1} $

I also will define the Fourier transform of a function depending on the variabl $p$ as
    $$ F^{-1}[f(p)](z) =\frac{1}{\sqrt{2\pi h}} \int_{-\infty}^{\infty} f(p) \exp\left\{\frac{ipz}{h}\right\} dp $$

For the moment I will only focus on the product of the Gaussian and the hermite polynomial with general $\alpha$ and $\beta$, and only later I will try to map the two parameters with the respective functions in the argument of theSTA wave function.

I want to compare three quantities:
- the analytical solution
- the numerical solution obtained by using the Fourier transform on the full STA wave function in monentum representation
- the numerical solution obtained by using the Fourier transform on only the spatial part of the STA wave function in monentum representation and then multiplying the time dependent one

In the following then I will assume that the STA wave function in momentum representation is of the form

$$ \frac{\sqrt{\xi_0}}{\sqrt[4]{\pi }
\exp \left(-i \left(n+\frac{1}{2}\right) \int_0^t \frac{\xi_0^2 U}{2 b(\tau )^2} \, d\tau \right) 
1/\sqrt{2^n n! b(t)}}
 e^{-p^2/2\alpha} \mathcal{H}_n(\beta p) $$

 with $\alpha$ and $\beta$ to be determined.
=#

"""
    `time_dependent(n,t, c::Control)`
    Return the solely time dependent part of the STA wave function in momentum representation.
"""
function time_dependent(n, t, c::Control)
    ξ0, U = EBJJ.scaling_ξ0(c), c.U # constants
    b(t) = auxiliary(t, c) # Auxiliary function
    normalisation = sqrt(EBJJ.scaling_ξ0(c)) / (π)^0.25 # Normalisation factor
    itp = EBJJ.interpolation_integral(c) # Interpolating function for the integral of the phase
    imaginary_phase(n, t) = exp(-im * (n + 0.5) * ξ0^2 * U / 2 * itp(t)) # Imaginary phase
    return normalisation * imaginary_phase(n, t) / sqrt(2^n * factorial(n) * b(t))
end
"""
    `spatial(n, p, α, β)`
Return the solely spatial part of the STA wave function in momentum representation in terms of the parameters α and β.
"""
function spatial(n, p, α, β)
    gaussian(p, α) = exp(-p^2 / (2 * α)) # Gaussian function
    hermite(n, p, β) = SpecialPolynomials.basis(Hermite, n)(β * p) # Hermite polynomial
    return gaussian(p, α) * hermite(n, p, β)
end
"""
    `analytical_solution(n, x, t, α, β)`
    Return the analytical solution of the Fourier transform of the product of a Gaussian function and a Hermite polynomial, where the argument of the two functions are two complex numbers `α` and `β`.
"""
function analytical_solution(n, x, h, α, β)
    γ = sqrt(Complex(2 * α * β^2 - 1)) # Auxiliary variable
    return im^n * sqrt(Complex(β)) * exp(-α * x^2 / 2) * γ^n * SpecialPolynomials.basis(Hermite, n)(α * β / γ * x)
end
"""
    `fourier_transform(x, h, f::Function)`
Give the Fourier transform of a function `f` depending on the variable `p` with constant factor `h`.
"""
function fourier_transform(x, h, f::Function)
    return quadgk(p -> f(p) * exp(-im * p * x / h), -10.0, 10.0, rtol=1e-3)[1] / sqrt(2 * π * h)
end

@testset "Fourier check, general α, β" begin
    α, β = 1.2, 0.1 # Parameters
    n = 2
    c = ControlFull()
    h = 1 / c.N
    t = [0.0, c.T / 2, c.T]
    full(p) = time_dependent(n, t[1], c) * spatial(n, p, α, β)
    analytic(x) = analytical_solution(n, x, h, α, β)
    xrange = range(-1, 1, length=2)
    full_x = [quadgk(p -> full(p) * exp(im * p * x), -10.0, 10.0, rtol=1e-3)[1] / sqrt(2 * π) for x in xrange]
    @test full_x ≈ analytic.(xrange)
    # @test analytic.(xrange) == 1.0
    # @test Spatial.(xrange) ≈ analytic.(xrange)
    # @test Full.(xrange) ≈ Spatial.(xrange)
end
