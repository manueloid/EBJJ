using QuadGK
using Interpolations
using BenchmarkTools
using ForwardDiff
using SpecialPolynomials
using EBJJ

#=
# Fourier checks
In this set of tests I am going to numerically check if the simplifications I made to turn the STA wave function from momentum to position representation are correct.

I have already checked in [the relevant repository](~/Repos/ExternalBJJ/fourier_transform.wlnb) that the Fourier transform of the product of a Gaussian function and of a Hermite polynomial has a specific close form.

Here I only want to check if the normalisation and the orthogonality hold.
This can be done by considering the fact that the Fourier transform of the whole STA wave function in momentum representation is the product of the Fourier transform of the spatial part and of the time dependent part.

Hence the STA wave function in position representation is just the product of the time dependent term and the Fourier transform of term dependent on the momentum varible $p$.

Plans:
First, I will define the analytical solution of the Fourier transform of the product of a Gaussian function of the form $ e^{-p^2/2\alpha}$ and a hermite polynomial of the form $ \mathcal_n(\beta p) $

It is given by:
$$ i^n \sqrt{\alpha} \exp\left\{\frac { -\alpha x^2}{2} \right\} \gamma^n \mathcal{H}_n\left( \frac{\alpha\beta}{\gamma}x\right)$$
where
$ \gamma = \sqrt{2 \alpha^2 \beta^2 - 1} $

Then I will define the whole STA wave function in position representation and I will check if the normalisation and the orthogonality hold.

I will start from the normalisation condition.
I will check the normalisation for two different kinds of functions, that in theory should give the same result:
1. the product of the STA wave function in position representation and its complex conjugate with no simplifations
2. the product of the STA wave function in position representation and its complex conjugate with the simplifications already implemented in the function definition.

For the moment I will just use the parameter $\alpha$ and $\beta$ to be complex numbers and I will check the normalisation for some fixed times.
If everything works according to plan, I will then check if the same relation holds if I substitute $\alpha$ and $\beta$ with the respective time depending functions
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
