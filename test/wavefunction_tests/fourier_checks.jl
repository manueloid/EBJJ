#=
# Fourier checks

In this file I am going to numerically check if the analytic form of the Fourier transform of the STA wave function is correct.
The plan is to numerically implement the Fourier transform and then pass the same value to the analytic form of the Fourier transform and check if the two are equal.

I will first do that for a general complex number and only later I will pass the actual value of the parameter $ \alpha^2$.

I already tested that the analytic solution can be expressed in terms of only on parameter $ \alpha^2$ and its real part $ R_\alpha$.
I will use this fact to simplify the code.
=#
#=
### 1 Fourier transform general complex number

The function I need to find the Fourier transform of is the following:
$$
    \exp\left\{	- \frac{p^2}{2} \alpha^2\right\}
    \mathcal{H}_n\left( p R_/alpha \right)
    \tag{1}
$$

where in this case the parameter $ \alpha $ is not dependent on time $ t $, this will be discussed later.
In the context of quantum mechanics, the Fourier transform is quite different as we need an extra term $ \hbar$ in the exponent and in the normalisation term.
In particular, the integral is
$$ 
	\frac{1}{\sqrt{2\pi\hbar}}
	\int_{-\infty}^{+\infty}
    \exp\left\{	- \frac{p^2}{2} \alpha^2\right\}
    \mathcal{H}_n\left( p R_/alpha \right)
	e^{\frac{i}{\hbar}px} dp
    \tag{2}
$$

And we can simplify this expression by defining a new variable  $ \tilde{z} = z / \hbar$.

The analytic solution of this integral is given by
$$
    \frac{i^n}{\sqrt{\hbar}} 
    \frac{1}{\alpha}
    \left[\frac{\alpha^{2*}}{\alpha^2}\right]^{n/2}
    \exp\left\{- \frac{\tilde{z}^2}{2\alpha^2}\right\}
    \mathcal{H}_n\left( \tilde{z} \frac{R_\alpha}{|\alpha^2} \right)
    \tag{3}
$$

I am going to test this relation for the new defined variable $ \tilde{z} $, I can always pass the actual value of $ z $ later.

#### 1.1 Code implementation

In terms of the implementation of the code, I will define two versions of each term of the function.
One will be the numerical solution to the integral and I will add the suffix `_num`, while the other one will be the analytical solution of the integral and I will use the suffix `_an`.

I will only care about the Gaussian term and the Hermite polynomial, as the other terms do not depend on the variable $ z $.
Moreover, I will only define the full function, I will not split the whole function in smaller chunks.
This is because this file is only a test and I do not need any optimisation.
=#

he(n::Int64, z::Float64) = SpecialPolynomials.basis(Hermite, n)(z)
"""
    toint(n::Int64, p::Float64, α2::ComplexF64)
This function is the integrand of the Fourier transform of the STA wave function.
It is basically only the product of the Gaussian term and the Hermite polynomial.
The function return the value of the integrand for the momentum `p` and for the parameter `α2`, given the excitation level `n`.
Is the one I am going to pass to the `quadgk` function.
"""
function toint(n::Int64, p::Float64, α2::ComplexF64)
    r = sqrt(real(α2))
    return exp(-p^2 * α2 / 2) * he(n, p * r)
end
"""
    fourier(n::Int64, z::Float64, α2::ComplexF64)
Numerically return the value of the Fourier transform of `toint` for the excitation level `n`, the position `z` and the parameter `α2`.
It also takes the value of the parameter `h` which is the Planck constant.
"""
function fourier_num(n::Int64, z::Float64, α2::ComplexF64, h::Float64)
    return quadgk(p -> toint(n, p, α2) * exp(im * p * z), -Inf, Inf, atol=1e-7)[1] / sqrt(2 * π * h)
end
"""
    fourier_an(n::Int64, z::Float64, α2::ComplexF64, h::Float64)
Return the value of the analytic solution of the Fourier transform of the STA wave function for the excitation level `n`, the position `z` and the parameter `α2`.
"""
function fourier_an(n::Int64, z::Float64, α2::ComplexF64, h::Float64)
    r = sqrt(real(α2))
    return im^n / sqrt(h) * 1 / sqrt(α2) * (conj(α2) / α2)^(n / 2) * exp(-z^2 / (2 * α2)) * he(n, z * r / abs(α2))
end
@testset "Fourier transform general complex" begin
    for _ in 1:10000
        n = rand(0:4)
        z = rand()
        α2 = rand(ComplexF64)
        h = 0.01
        @test isapprox(fourier_num(n, z, α2, h), fourier_an(n, z, α2, h), atol=1e-4)
    end
end


normalisation(η::ComplexF64, n::Int64) = (real(η) / pi)^(1 / 4) / sqrt(2^n * factorial(n)) * im^n * sqrt(conj(η) / η)^n / sqrt(η)
gaussian(z::Float64, η::ComplexF64) = exp(-z^2 / (2 * η))
hermite(n::Int64, z::Float64, η::ComplexF64) = SpecialPolynomials.basis(Hermite, n)(z * sqrt(real(η)) / abs(η))
spatial_fourier(n::Int64, z::Float64, η::ComplexF64) = normalisation(η, n) * gaussian(z, η) * hermite(n, z, η)

@testset "normalisation STA position general" begin
    η = rand(ComplexF64)
    k = rand(Float64)
    f(k, n) = exp(im * (n + 1 / 2) * k)
    for n in 0:5
        @test isapprox(quadgk(x -> conj(f(k, n) * spatial_fourier(n, x, η)) * f(k, n) * spatial_fourier(n, x, η), -Inf, Inf, atol=1e-7)[1], 1.0, atol=1e-3)
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
    α(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) - im * db(t) / (U * b(t)) # Parameter of the Gaussian term
    imag_phase_integrand(t::Float64) = sqrt(2J0 * N * U) / b(t)^2
    φ(t::Float64) = quadgk(τ -> imag_phase_integrand(τ), 0.0, t)[1]
    return exp(-im * (n + 1 / 2) * φ(t)) * spatial_fourier(n, x, α(t))
end
@testset "normalisation with system parameters" begin
    c = ControlFull(10, 0.03)
    norm(n, t) = quadgk(x -> conj(wave_function(n, t, x, c)) * wave_function(n, t, x, c), -Inf, Inf, atol=1e-7)[1] |> real
    for t in range(0.0, c.T, length=100), n in 0:5
        @test isapprox(norm(n, t), 1.0, atol=1e-3)
    end
end
@testset "Orthogonalisation with system parameters" begin
    c = ControlFull(10, 0.03)
    orth(n, m, t) = quadgk(x -> conj(wave_function(n, t, x, c)) * wave_function(m, t, x, c), -Inf, Inf, atol=1e-7)[1] |> real
    for _ in 1:1000
        t, n, m = rand(0.0:c.T), rand(0:3), rand(4:7)
        if n == m
            @test isapprox(orth(n, m, t), 1.0, atol=1e-3)
        else
            @test isapprox(orth(n, m, t), 0.0, atol=1e-3)
        end
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
    α(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) + im * db(t) / (U * b(t)) # Parameter of the Gaussian term
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
function test_ground_state(η::ComplexF64, n::Int64)
    lhs(x) = lhs_spatial_fourier(n, x, conj(η))
    rhs(x) = gaussian(x, η)
    return quadgk(x -> lhs(x) * rhs(x), -Inf, Inf, atol=1e-7)[1]
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
    α(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) - im * db(t) / (U * b(t)) # Parameter of the Gaussian term
    αc(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) + im * db(t) / (U * b(t)) # Complex Conjugate of the Parameter of the Gaussian term 
    lhs(x) = lhs_spatial_fourier(n, x, αc(t))
    rhs(x) = gaussian(x, α(t))
    return quadgk(x -> lhs(x) * rhs(x), -Inf, Inf, atol=1e-7)[1]
end

@testset "Testing lhs actual values" begin
    c = ControlFull(10, 0.3)
    for t in range(0.0, c.T, length=100), n in 0:5
        if n == 0
            @test isapprox(test_ground_state(n, t, c), 1.0, atol=1e-3)
        else
            @test isapprox(test_ground_state(n, t, c), 0.0, atol=1e-3)
        end
    end
end
