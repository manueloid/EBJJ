#=
# Testing the validity of the single parameter

In this file I am going to check if all the calculations I made in the derivation of the STA wave function with only on complex parameter, is valid.

The plan is to define the wave function with all the different paremeters and then check if it is equal to the one where only one parameter is passed.

The following is the wave function with only one parameter.
=#

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
