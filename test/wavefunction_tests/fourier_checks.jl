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

I also will define the Fourier transform of a function depending on the variable $p$ as
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

For the moment I will only focus on the spatial term of the product of the Gaussian term and the Hermite polynomial, I will then check if the normalisation and the orthogonality hold and finally I will define the analytic solution of the Fourier transform of the product and test the results.
Again, the Gaussian and the Hermite polynomial will be in terms of the two parameters $\alpha$ and $\beta$.
=#

"""
    `spatial(n::Int64, p::Float64, α::ComplexF64, β::Float64)`
Return the product of the Gaussian function and the Hermite polynomial in momentum representation and in terms of the parameters `α` and `β`.
"""
function spatial(n::Int64, p::Float64, α::ComplexF64, β::Float64)
    return exp(-p^2 / (2 * α^2)) * SpecialPolynomials.basis(Hermite, n)(β * p)
end
"""
    `analytic(n::Int64, x::Float64, α::ComplexF64, β::Float64)`
Return the analytic solution of the Fourier transform of the product of the Gaussian function and the Hermite polynomial in position representation and in terms of the parameters `α` and `β`.
"""
analytic(n) = n

@testset "Testing normalisation momentum for α and β" begin
    α = 1.0 + 0.0 * im
    β = 1.0
    @test quadgk(
        x -> abs(
            spatial(2, x, 1.0 + 0.0 * im, 1.0) *
            conj(spatial(0, x, 1.0 + 0.0 * im, 1.0))
        ), -Inf, Inf)[1] ≈ sqrt(π)
end














































