using Test
#=
# Testing the simplifications

In this file I am goint to see if the simplifications I made when calculating the quantities of the form
$$
 \langle \chi_n(z,t) | \hat{O} | \chi_0(z,t) \rangle
 \tag{1}
$$
where $ \hat{O} $ is a generic operator, are correct.

In particular, I know that the operators acts only on the spatial part of the wave function, so I can simplify the expressions by considering all the terms only depending on the time variable $ t $.

I have already tested that:
1- The STA wave function can be parametrized using only one complex value $ \alpha^2 $ and its real part $ R_\alpha $.
2- The analytic solution of the STA Fourier transform is correct.

So in the following I am going to use these two facts to simplify the expressions of the quantities of the form (1).
=#
#=
## 1 Whole function

The whole STA wave function, result of the Fourier transform of the solution of the Harmonic oscillator in momentum representation, is given by:
$$
    \chi_n(z,t) =
    \left[\frac{R_\alpha^2}{\pi}\right]^{1/4}
    \frac{1}{\sqrt{2^n n!}} % Hermite polynomial normalisation factor   
	\exp\left\{-i(n + 1/2)\ \int_{0}^{t}d \tau ~ U R_\alpha^2(\tau)\right\} % time dependent phase factor
    \frac{i^n}{\sqrt{\hbar}} \frac{1}{\alpha} \left(\frac{\alpha^{2*}}{\alpha^2}\right)^{n/2} % normalisation factor
    \exp\left\{-\frac{z^2}{2 \hbar^2 \alpha^2} \right\} % Gaussian term
    \mathcal{H}_n\left(\frac{R_\alpha}{|\alpha^2|}~\frac{z}{\hbar}\right) 
    \tag{2}
$$
where $ \alpha $ is a complex number and $ R_\alpha $ is its real part.

I am going to define the following functions to make the code more readable.
All of them will have the suffix `_wh` to indicate that they are the whole function.

- `norm_wh(n::Int64, α2::ComplexF64, h::Float64)`: Normalisation of the whole wave function
This corresponds to the following formula:
$$ 
    \left[\frac{R_\alpha^2}{\pi}\right]^{1/4}
   \frac{1}{\sqrt{2^n n!}} % Hermite polynomial normalisation factor   
    \frac{i^n}{\sqrt{\hbar}} \frac{1}{\alpha} \left(\frac{\alpha^{2*}}{\alpha^2}\right)^{n/2} % normalisation factor
    \tag{3}
$$

- `ga_wh(z::Float64, α2::ComplexF64)`: Gaussian term of the whole wave function
This corresponds to the following formula:
$$
    \exp\left\{-\frac{z^2}{2 \hbar^2 \alpha^2} \right\} % Gaussian term
    \tag{4}
$$

-  `he_wh(n::Int64, z::Float64, α2::ComplexF64)`: Hermite polynomial of the whole wave function
This corresponds to the following formula:
$$
    \mathcal{H}_n\left(\frac{R_\alpha}{|\alpha^2|}~\frac{z}{\hbar}\right) 
    \tag{5}
$$

- `ph_wh` : Phase factor of the whole wave function
This corresponds to the following formula:
$$
    \exp\left\{-i(n + 1/2)\ \int_{0}^{t}d \tau ~ U R_\alpha^2(\tau)\right\} % time dependent phase factor
    \tag{6}
$$

I am going to implement `ph_wh` only later as I need to test the normalisation first for a general complex number and here I do not explicitely define the time dependent factor $ \R_\alpha $
=#
#=
##  2 Normalisation

I need to test the normalisation so I need the function $ |\chi_n(z,t)|^2 $ and then I am going to evaluate the integral   $ \int_{-\infty}^{+\infty} |\chi_n(z,t)|^2 dz $.

For the moment I am not going to perform many simplifications, I am going to define the whole function and integrate.
The only simplification I am going to make is the one relative to the complex phase integral as I assume that it will vanish when pairing the wave function with its complex conjugate.

As usual, I will pass a general complex number and only later I will pass the actual values of the system.
=#

he(n::Int64, z::Float64) = SpecialPolynomials.basis(Hermite, n)(z) # one liner to speed up the writing of the polynomial
"""
    norm_wh(n::Int64, z::Float64, α2::ComplexF64)
Return the normalisation of the whole wave function for excitation `n` and for a general complex number `α2`.
"""
function norm_wh(n::Int64, α2::ComplexF64, h::Float64)
    r = sqrt(real(α2))
    return (r^2 / pi)^(1 / 4) / sqrt(2^n * factorial(n)) * im^n / sqrt(h) * 1 / sqrt(α2) * (conj(α2) / α2)^(n / 2)
end
"""
    ga_wh(z::Float64, α2::ComplexF64)
Return the Gaussian part of the wave function at a given position for a general complex number `α2`.
The idea is to pass the variable `z/h` when calling the function, so there is no need to pass the value `h` in this functio
"""
function ga_wh(z::Float64, α2::ComplexF64)
    return exp(-z^2 / (2 * α2))
end
"""
    he_wh(n::Int64, z::Float64, α2::ComplexF64)
Return the valuee of the Hermite polynomial for excitation number `n` at the position `z` for a general complex number `α2`.
Again, no need to pass the value `h` here, pass the value `z/h` when calling the function.
"""
function he_wh(n::Int64, z::Float64, α2::ComplexF64)
    return he(n, z * sqrt(real(α2)) / abs(α2))
end
"""
    whole(n::Int64, z::Float64, α2::ComplexF64, h::Float64)
Return the value of the whole wave function for excitation number `n` at the position `z` for a general complex number `α2`.
In this case, the value `h` is needed and passed internally to the relevant functions.
"""
function whole(n::Int64, z::Float64, α2::ComplexF64, h::Float64)
    return norm_wh(n, α2, h) * ga_wh(z / h, α2) * he_wh(n, z / h, α2)
end
@testset "Normalisation test general complex number" begin
    n = rand(0:4)
    h = 0.1
    for _ in 1:10000
        α2 = rand(ComplexF64)
        normalisation(z::Float64) = abs(whole(n, z, α2, h))^2
        res = quadgk(normalisation, -Inf, Inf,)[1]
        @test isapprox(res, 1.0, atol=1e-4)
    end
end



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
    α(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) - im * db(t) / (U * b(t)) # Parameter of the Gaussian term
    imag_phase_integrand(t::Float64) = sqrt(2J0 * N * U) / b(t)^2
    φ(t::Float64) = quadgk(τ -> imag_phase_integrand(τ), 0.0, t)[1]
    return exp(-im * (n + 1 / 2) * φ(t)) * spatial_fourier_wh(n, x, α(t))
end
ground_state(t, x, c::Control) = wave_function(0, t, x, c)
pairing(n::Int64, t, x, c::Control) = ground_state(t, x, c) * conj(wave_function(n, t, x, c))
normalisation(η::ComplexF64, n::Int64) = (real(η) / pi)^(1 / 2) / sqrt(2^n * factorial(n)) * (-im)^n * sqrt(conj(η) / η)^n / abs(η)
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

@testset "Testing the simplifications" begin
    c = ControlFull()
    for _ in 1:1000, n in 0:5
        t, z = rand(Float64), rand(Float64)
        @test pairing(n, t, z, c) ≈ simplified(n, t, z, c)
    end
end

#=
### Full normalisation test

If the normalisation holds, the integral over the whole time interval should yield the width of the time interval.
I am going to test this.
=#

@testset "Normalisation over time interval" begin
    tf = rand()
    np = rand(10:10:100)
    c = ControlFull(np, tf)
    lim = 1e3
    for n in 2:2:4
        wf(x) = conj(wave_function(n, x[1], x[2], c)) * wave_function(n, x[1], x[2], c) |> real
        result = hcubature(wf, [0.0, -lim], [tf, lim], rtol=1e-11, atol=1e-10)[1]
        @test result ≈ tf
    end
end
