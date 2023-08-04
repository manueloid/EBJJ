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
If everything works according to plan, I will then check if the same relation holds if I substitute $\alpha$ and $\beta$ with the respective time depending functions.

Since calculating the normalisation in the case of general $\alpha$ and $\beta$ is quite cumbersome, I will only check the orthogonality. If that holds, I will then move on to define $\alpha$ and $\beta$ as functions of time.
=#

"""
    `analytic(n::Int64, x::Float64, α::ComplexF64, β::Float64)`
Return the value of the analytic solution of the Fourier transform of the product of a Gaussian function and a Hermite polynomial, at the point `x`.
"""
function analytic(n::Int64, x::Float64, α::ComplexF64, β::Float64)
    γ = sqrt(2.0 * α^2 * β^2 - 1.0)
    return im^n * α * exp(-α^2 * x^2 / 2.0) * γ^n * SpecialPolynomials.basis(Hermite, n)(x * α^2 * β / γ)
end
"""
    `orthogonality_check(n::Int64, m::Int64, α::ComplexF64, β::Float64)`
Take the integral over the variable `x` of two analytic solutions of the Fourier transform of the product of a Gaussian function and a Hermite polynomial, with different values of `n` and `m`.
"""
function orthogonality_check(n::Int64, m::Int64, α::ComplexF64, β::Float64)
    return quadgk(x -> conj(analytic(n, x, α, β)) * analytic(m, x, α, β), -1e9, 1e9)[1]
end
@testset "Orthogonality test, general α and β" begin
    α = 0.4 + 0.3im
    β = 0.07
    # Check that the integral is zero if n != m
    for n in 0:5, m in 0:5
        if n != m
            @test orthogonality_check(n, m, α, β) ≈ 0.0 + 0.0im
        end
    end
end

#=
The test have been passed, so now I can use the actual values of $\alpha$ and $\beta$ and not some general ones.
In particular, we have that:
1. $$ \alpha^2 = \left(- \frac{\xi_0^2}{b} + \frac{2i\dot{b}}{Ub} \right)^-1 $$
2. $$ \beta = \frac{\xi_0}{b} $$

The plan is to test both the orthogonality and the normalisation.
I think the best plan now is to define two functions that will return the values of $\alpha $ and $\beta$ at a given time and then pass it to the function `analytic`.

First I will check the orthogonality as it easier to do so, if that goes well I will continue with the normalisation, and implement multiple checks for it.

I think the best way to implement the code is to put everything into a big function, as I need to define the auxiliary function $b(t)$, and it is not ideal to have the code redefine it every time.
I will check if there is any more flexible and performant way to perform this task, but for the moment I only want to get the job done.
=#
"""
    `wave_function(n, t, x, c::Control)`
Return the value of the STA wave function in position representation at the point `x` and at time `t`, given the control parameters `c`.
This is the function with absolutely no simplifications and it is the one I am going to test against.
"""
function wave_function(n, t, x, c::Control)
    ξ0, U = EBJJ.scaling_ξ0(c), c.U # constants
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    normalisation = sqrt(EBJJ.scaling_ξ0(c)) / (π)^0.25 # Normalisation factor
    itp = EBJJ.interpolation_integral(c) # Interpolating function for the integral of the phase
    imaginary_phase(n, t) = exp(-im * (n + 0.5) * ξ0^2 * U / 2 * itp(t)) # Imaginary phase
    α2(t) = (-ξ0^2 / b(t)^2 + 2im * db(t) / (U * b(t)))-1 # α^2 term to be used in the gaussian
    β(t) = ξ0 / b(t) # β term to be used in the hermite polynomial
    spatial(n,t,x) = analytic(n, x, sqrt(α2(t)), β(t)) # Spatial term
    excitation(n, t) = sqrt(2^n * factorial(n) * b(t)) # Excitation term
    return normalisation * excitation(n, t) * imaginary_phase(n, t) * spatial(n, t, x)
end
"""
