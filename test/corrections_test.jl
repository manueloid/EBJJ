using Test
#=
# Testing the corrections

In this test suite, the plan is to set up two different ways of evaluating the corrections, to numerically check if the simplifications cause any disruption.

I will first check if the two methods give the same results without calculating the corrections, then I will implement the code to calculate the corrections to see if the two methods are equivalent.

I am going to reuse the definition of the functions I used in different tests.
There is no need to use a general complex value $ \eta $, as I already checked the correctness of the code implementation.
I will thus start by using the function definition that takes into account the parameters of the system.
=#

"""
    argument(t, c::Control)
Return the argument of the Gaussian part of the STA wave function at a time `t` given the control parameter `c`.
"""
function argument(t, c::Control)
    J0, N, U = c.J0, c.N, c.U
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    α(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) - im * db(t) / (U * b(t)) # Parameter of the Gaussian term
    return α(t)
end
normalisation_wh(n::Int64, η::ComplexF64) = (real(η) / pi)^(1 / 4) / sqrt(2^n * factorial(n)) * im^n * sqrt(conj(η) / η)^n / sqrt(η)
normalisation_wh(n::Int64, t::Float64, c::Control) = normalisation_wh(n, argument(t, c))

gaussian_wh(z::Float64, η::ComplexF64) = exp(-z^2 / (2 * η))
hermite_wh(n::Int64, z::Float64, η::ComplexF64) = SpecialPolynomials.basis(Hermite, n)(z * sqrt(real(η)) / abs(η))
spatial_fourier_wh(n::Int64, z::Float64, η::ComplexF64) = normalisation(η, n) * gaussian(z, η) * hermite(n, z, η)
"""
     wave_function(n::Int64, t, x, c::Control)
Return a general STA wave function for excitation `n` at time `t` and position `z`.
This function is helpful to debug code but is not the best when it comes to performances
"""
function wave_function(n::Int64, t, x, c::Control)
    J0, N, U = c.J0, c.N, c.U
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    α(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) - im * db(t) / (U * b(t)) # Parameter of the Gaussian term
    imag_phase_integrand(t::Float64) = sqrt(2J0 * N * U) / b(t)^2
    φ(t::Float64) = quadgk(τ -> imag_phase_integrand(τ), 0.0, t)[1]
    return exp(-im * (n + 1 / 2) * φ(t)) * spatial_fourier_wh(n, x, α(t))
end

ground_state(t, x, c::Control) = wave_function(0, t, x, c)
pairing(n::Int64, t, x, c::Control) = ground_state(t, x, c) * conj(wave_function(n, t, x, c))
normalisation(n::Int64, η::ComplexF64) = (real(η) / pi)^(1 / 2) / sqrt(2^n * factorial(n)) * (-im)^n * sqrt(conj(η) / η)^n / abs(η)
normalisation(n::Int64, t, c::Control) = normalisation(n, argument(t, c))

gaussian(z::Float64, η::ComplexF64) = exp(-z^2 / (2 * η))
hermite(n::Int64, z::Float64, η::ComplexF64) = SpecialPolynomials.basis(Hermite, n)(z * sqrt(real(η)) / abs(η))
spatial_fourier(n::Int64, z::Float64, η::ComplexF64) = normalisation(η, n) * gaussian(z, η) * hermite(n, z, η)
"""
    simplified(n::Int64, t, x, c::Control)
Return the value of the product between complex conjugate of a STA wave function and the ground state, assuming all the simplifications have been carried out.
"""
function simplified(n::Int64, t, x, c::Control)
    J0, N, U = c.J0, c.N, c.U
    # Definition of the auxiliary functions
    b(t) = auxiliary(t, c)
    db(t) = ForwardDiff.derivative(b, t)
    α(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) - im * db(t) / (U * b(t))  # Parameter of the Gaussian term
    αc(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) + im * db(t) / (U * b(t)) # Complex Conjugate of the Parameter of the Gaussian term 
    imag_phase_integrand(t::Float64) = sqrt(2J0 * N * U) / b(t)^2
    φ(t::Float64) = quadgk(τ -> imag_phase_integrand(τ), 0.0, t)[1]
    return exp(im * n * φ(t)) * spatial_fourier(n, x, αc(t)) * gaussian(x, α(t))
end

@testset "Testing the normalisation" begin
    for _ in 1:1000
        n, tf = rand(0:2:6), rand(Float64)
        c = ControlFull(n, tf )
        t, z = rand(Float64) * c.T, rand(Float64)
        @test isapprox( conj(normalisation_wh(n, t, c)) * normalisation_wh(0, t, c), normalisation(n, t, c), atol = 1e-4) 
    end
end
