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


### 1 Fourier transform general complex number
First, I will define the analytical solution of the Fourier transform of the product of a Gaussian function of the form $ e^{-p^2/2\alpha}$ and a hermite polynomial of the form $ \mathcal_n(\beta p) $

It is given by:
$$ i^n \sqrt{\alpha} \exp\left\{\frac { -\alpha x^2}{2} \right\} \gamma^n \mathcal{H}_n\left( \frac{\alpha\beta}{\gamma}x\right)$$
where
$ \gamma = \sqrt{2 \alpha^2 \beta^2 - 1} $

I am going to use the fact that if we call the argument of the Gaussian part $\eta^2$ (such that $ 1/\alpha^2 = \eta^2$), then the parameter $ \beta $ is nothing more than $ \sqrt{\Re{\eta^2}} $.
Thus, we can simplify the term $ \alpha ^2 \beta ^2 -1 $ that appears twice in the Fourier transform of the product of the Gaussian term and the Hermite polynomial as $ \eta^{2*} / \eta^2 $.
With this simplification, the term $\gamma$ becomes $ \eta^{2*} / \eta^{2} $.
Moreover, the term $\alpha^2 \beta $ appearing in the numerator of the argument of the Hermite polynomial becomes $ \sqrt{\Re{\eta^2}} \eta^{2} $.

This simplification is helpful as I can now write the Fourier transform of the product of the Gaussian term and the Hermite polynomial in term of only one parameter, $\eta$, which is the argument of the Gaussian term of the STA wave function in momentum representation.

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
    `analytic(n::Int64, x::Float64, η::ComplexF64)`
Return the value of the analytic solution of the Fourier transform of the product of a Gaussian function and a Hermite polynomial, at the point `x`, given the complex parameter `η` and the order of the Hermite polynomial `n`.
If we want to connect this to the STA wave function, the parameter `η` is actually `η²`.
"""
function analytic(n::Int64, x::Float64, η::ComplexF64)
    γ = sqrt(conj(η) / η) # this is the √α²β² - 1 factor
    num = sqrt(real(η)) / η # This is the numerator inside the Hermite polynomial
    return im^n / sqrt(η) * γ^n * exp(-x^2 / (2 * η)) * SpecialPolynomials.basis(Hermite, n)(num * x / γ)
end
"""
    `orthogonality_check(n::Int64, m::Int64, α::ComplexF64, β::Float64)`
Take the integral over the variable `x` of two analytic solutions of the Fourier transform of the product of a Gaussian function and a Hermite polynomial, with different values of `n` and `m`, given a complex parameter `η`.
"""
function orthogonality_check(n::Int64, m::Int64, η::ComplexF64)
    return quadgk(x -> conj(analytic(n, x, η)) * analytic(m, x, η), -1e8, 1e8, atol=1e-3)[1]
end
@testset "Orthogonality test, general complex number η" begin
    η = 0.4 + 0.2im
    # Check that the integral is zero if n != m
    for n in 0:1, m in 0:2
        if n != m
            @test orthogonality_check(n, m, η) ≈ 0.0 + 0.0im
        end
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
