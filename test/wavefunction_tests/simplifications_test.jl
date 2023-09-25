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
=#
#=
### 2.1 Normalisation with general complex parameter
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
    h = 0.1
    for _ in 1:10000
        n = rand(0:4)
        α2 = rand(ComplexF64)
        normalisation(z::Float64) = abs(whole(n, z, α2, h))^2
        res = quadgk(normalisation, -Inf, Inf,)[1]
        @test isapprox(res, 1.0, atol=1e-4)
    end
end

#=
### 2.2 Normalisation with real parameter

In this section I am going to pass the actual parameters of the system and check if the normalisation is conserved.
I will define a function that takes the time as a parameter an returns the value of $ \alpha^2 $ at that time.
I will then pass said value to the functions defined above.

I will first try to solve this for different times and then I will integrate all over the time interval to see if the value of the integral is equal to the width of the interval 
=#

function α2(t::Float64, c::Control)
    J0, N, U = c.J0, c.N, c.U
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    α(t::Float64) = (sqrt(8J0 * N / U) / b(t)^2 - 2im / U * db(t) / b(t))
    return α(t)
end
@testset "Normalisation test system parameters" begin
    tf = rand()
    np = rand(10:10:50)
    c = ControlFull(np, tf)
    h = 1 / np
    for _ in 1:10000
        n = rand(0:4)
        t = rand(0.0:tf)
        normalisation(z::Float64) = abs(whole(n, z, α2(t, c), h))^2
        res = quadgk(normalisation, -Inf, Inf,)[1]
        @test isapprox(res, 1.0, atol=1e-4)
    end
end

@testset "Normalisation system whole interval" begin
    tf = rand()
    np = rand(10:10:50)
    c = ControlFull(np, tf)
    h = 1 / np
    for n in 0:4
        normalisation(v) = abs(whole(n, v[1], α2(v[2], c), h))^2
        res = hcubature(normalisation, [-10.0, 0.0], [10.0, tf])[1]
        @test isapprox(res, tf, atol=1e-4)
    end
end
