#=
# Single Parameter testing

In the notes I already set up the calculations where only one parameter is used.

Just for reference, the formula of the Fourier transform of the STA wave function where the argument of the Gaussian part is set as the parameter $ \alpha $ is given by
$$
    \chi_n(z,t) = \left [ \frac{R_\alpha^2}{π}]^{1/4} \frac{1}{\sqrt{2^nn!} \exp\left\{- i \left(n + \frac{1}{2} \right) \int_0^t d\tau ~ U R_\alpha^2(\tau) \right\}
    \frac{i^n}{\sqrt{\hbar}}\frac{1}{\alpha} \left[\frac{(\alpha^2)^*}{\alpha^2}\right]^{n/2}
    \exp{-\frac{z^2}{2 \hbar^2 \alpha^2}}
    \mathcal{H}_n\left(\frac{R_\alpha}{|\alpha^2|} \frac{z}{\hbar}\right)
    \tag{1}
$$

where - in terms of the external bosonic Josephson Junction parameters - the value $ \alpha ^2 $ is defined as
$$
  \alpha^2 = \frac{ \sqrt{8 J_0 N / U}}{b^2} - \frac{2 i \dot{b}}{U b} 
  \tag{2}
$$

and the term $ R $ is defined as the square root of the real part of $ \alpha ^2 $, i.e.
$$
R_\alpha = \sqrt{\Re(\alpha^2)} = \left(\frac{8J_0 N}{U}\right)^{1/4} \frac{1}{b}
  \tag{3}
$$

What I want to test here is if the simplifications I made when dealing with the complex numbers are correct.
In particular, I first want to check if the following equality holds
$$
2 \frac{R_\alpha^2}{\alpha^2} - 1 = /frac{(\alpha^2)^*}{\alpha^2} 
\tag{4}
$$

Then I want to numerially check if the following relation is true
$$
    \frac{R_\alpha}{\alpha^2 \sqrt{2 \frac{R_\alpha^2}{\alpha^2} - 1}} = \frac{R_\alpha}{|\alpha^2|}
\tag{5}
$$

As I said I am going to do that numerically, first I will pass some random complex number to see if the relations are correct and only then I will pass the actual values of $ \alpha^2 $ for different times $t$.
=#
#=
## 1 - Testing the relations with random complex numbers

Here I am going to define some functions to test the relations I discussed earlier, but in this case the argument will be only a general complex number.

The plan is to test the relation for a certain number of general complex number and make sure that the relation holds.
=#

using Test
"""
    test_1(param::ComplexF64)
Given a complex number `param`, this function tests if the relation (4) holds.
The relation is 2 ℜ(α) / α - 1 = (α*) / α
"""
function test_1(param::ComplexF64)
    r = sqrt(real(param))
    result = 2 * r^2 / param - 1 - (conj(param) / param)
    return result
end
@testset "testing relation 1" begin
    for _ in 1:10000
        @test isapprox(test_1(), 0.0, atol=1e-10)
    end
end
"""
    test_2(param::ComplexF64)
Given a complex number `param`, this function tests if the relation (5) holds.
The relation is α^2 * sqrt(2 * R_α^2 / α^2 - 1) = |α^2|
"""
function test_2(param::ComplexF64)
    result = (param * sqrt(conj(param)/param)) - abs(param)
    return result
end
@testset "testing relation 2" begin
    for _ in 1:10000
        param = 1000*rand(ComplexF64)
        @test isapprox(test_2(param), 0.0, atol=1e-10)
    end
end

#=
## 2 - Testing the relations with the actual values of α^2

Now it is time to pass the actual values of $ \alpha^2 $ and see if the relations hold.
There is no need to define any new functions, but only to get the value of $ \alpha^2 $ at a certain time $ t $ and then pass it to the functions I defined earlier.

The function α2 is the one that returns the value of $ \alpha^2 $ at a certain time $ t $.
This function will also accept some control parameter variable that are the ones that give the parameters of the system.

I would like to point out that in this case I did not make the substitution $ \hbar \rigtarrow 1 / N$.
=#

function α2(t::Float64, c::Control)
    J0, N, U = c.J0, c.N, c.U
    h = 1 / N
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    α(t::Float64) = ( sqrt(8J0 / (U * N)) / b(t)^2 - 2im / (U * N) * db(t) / b(t) ) / h 
    return α(t)
end
c = ControlFull()
t = range(0.0, c.T, length=10000)
@testset "testing relation 1" begin
    for i in t
        @test isapprox(test_1(α2(i, c)), 0.0, atol=1e-10)
    end
end
@testset "testing relation 2" begin
    for i in t
        @test isapprox(test_2(α2(i, c)), 0.0, atol=1e-10)
    end
end







normalisation_wh(η::ComplexF64, n::Int64) = (real(η) / pi)^(1 / 4) / sqrt(2^n * factorial(n)) * im^n * sqrt(conj(η) / η)^n / sqrt(η)
gaussian_wh(z::Float64, η::ComplexF64) = exp(-z^2 / (2 * η))
hermite_wh(n::Int64, z::Float64, η::ComplexF64) = SpecialPolynomials.basis(Hermite, n)(z * sqrt(real(η)) / abs(η))
spatial_fourier_wh(n::Int64, z::Float64, η::ComplexF64) = normalisation_wh(η, n) * gaussian_wh(z, η) * hermite_wh(n, z, η)
"""
     wave_function(n::Int64, t, x, c::Control)
Return a general STA wave function for excitation `n` at time `t` and position `z`.
This function is helpful to debug code but is not the best when it comes to performances
"""
function wave_function(n::Int64, t, x, c::Control)
    J0, N, U = c.J0, c.N, c.U
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    α(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) - 2*im * db(t) / (U * b(t)) # Parameter of the Gaussian term
    imag_phase_integrand(t::Float64) = sqrt(2J0 * N * U) / b(t)^2
    φ(t::Float64) = quadgk(τ -> imag_phase_integrand(τ), 0.0, t)[1]
    return exp(-im * (n + 1 / 2) * φ(t)) * spatial_fourier_wh(n, x, α(t))
end
ground_state(t, x, c::Control) = wave_function(0, t, x, c)

#=
Now I am going to define the wave function with all the parameters and then check if it is equal to the one with only one parameter.

I need to use the solution that I found on the internet about the Fourier transform of the product of a Gaussian and a Hermite polynomial.
This is given by the following equation:

$$
\begin{align}
	\exp\left\{-\frac{k^{2}}{2\alpha^{2}}\right\} \mathcal{H}_{n}(\betak)  e^{ikx} = \nonumber \\ % Integrand expression
	i^{n}\alpha \left(2\beta^{2}\alpha^{2} - 1\right)^{n / 2} % Factor in front of the Hermite polynomial
	\exp\left\{-\frac{x^{2}\alpha^{2}}{2}\right\} % Gaussian term
	\mathcal{H}_{n}\left(\frac{\alpha^{2}\betax}{\sqrt{2\alpha^{2}\beta^{2} - 1}}\right) ~. % Hermite polynomial 
	\tag{1}
\end{align}
$$

In this case, we have that $ a^2 = \left(2 \frac{\sqrt{2J_0 N / U}}{b^2}  - 2 \frac{i \dot{b}}{Ub}\right)^{-1}$ and $ \beta = \left(\frac{8 J_0 N}{U}\right)^{1/4}\frac{1}{b}
and
=#

function wave_function_full(n::Int64, t, z, c::Control)
    J0, N, U = c.J0, c.N, c.U
    k = 8J0 * N / U
    ω0 = sqrt(2J0 * N * U)
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    α(t::Float64) = ( sqrt(k) / b(t)^2 - 2im / U * db(t) / b(t) )^(-1) # Parameter of the Gaussian term
    β(t::Float64) = ( k )^(1/4) / b(t)
    imag_phase_integrand(t::Float64) = ω0 / b(t)^2
    φ(t::Float64) = quadgk(τ -> imag_phase_integrand(τ), 0.0, t)[1]
    # Definition of the whole function
    result(t,z) = 
    (sqrt(k) / pi)^(1 / 4) *                            # Normalisation constant factor 
    exp(-im * (n + 1 / 2) * φ(t)) *                     # Imaginary phase
    (2^n * factorial(n) * b(t))^(-1 / 2) *              # Normalisation factor depending on the excitation
    im^n * (2 * α(t)^2 * β(t)^2 - 1)^(n/2) * sqrt(α(t))            # Factor in front of the Hermite polynomial
    exp(-z^2 * α(t)^2 / 2) *                            # Gaussian term
    SpecialPolynomials.basis(Hermite, n)(α(t)^2 * β(t) * z / sqrt(2 * α(t)^2 * β(t)^2 - 1)) # Hermite polynomial
    return result(t,z)
end

norm_full(n::Int64, t, c::Control) = quadgk(z -> abs(wave_function_full(n, t, z, c))^2, -10000, 10000)[1]
norm(n::Int64, t, c::Control) = quadgk(z -> abs(wave_function(n, t, z, c))^2, -Inf, Inf)[1]
c = ControlFull()
t = 0.0
n = rand(0:10)
println("The norm of the wave function with only one parameter is $(norm_full(n, t, c))")
