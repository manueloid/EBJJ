#=
# Testing Schrodinger equation

I want to see if the Fourier transform of the STA wave function is a solution of the Schrodinger equation in position representation.

The Hamiltonian in position representation is given by
$$
    H = - 2 h^2 J(t) \partial_z^2 + \frac{UN}{4}z^2
    \tag{1}
$$
where $ J(t) $ is the control function and $ h $ is somethign resembling the Planck constant.
In this case this value is $ 1/N$, where $ N $ is the number of particles.

The Schrodinger equation can thus be written as
$$
    i \partial_t \chi(z,t) = H \chi(z,t)
    \tag{2}
$$
where $ \chi(z,t) $ is the STA wave function in position representation and already derived elsewhere.

I am going to use automatic differentiation, hopefully it will work.

The problem as always is the time dependent phase, that I cannot differentiate.
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


