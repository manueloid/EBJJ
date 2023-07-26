#=
# STA Wave function 

I need to things slowly and carefully, and write tests along the way.
In this file I will only focus on defining the wave function, in a different one I will actually evaluate the eSTA corrections.

## 1. Constant factor and helpful functions

First I will start by defining the constant factor $\xi_0$ and the functions $\alpha(t)$ and $\beta(t)$ that go into the Fourier transform of the STA wave function.
The names of the function will be (respectively) `scaling_ξ0`, `gaussian_arg` and `hermite_arg`.
=#

scaling_ξ0(c::Control) = sqrt(8c.J0 * c.N / c.U)
"""
    `gaussian_arg(t, c::Control)` 
    This function maps the argument of the Gaussian exponential of the STA wave function to the argument of the standard Gaussian function used to evaluate the Fourier transform of the product between a Gaussian function and a Hermite polynomial.
"""
function gaussian_arg(t, c::Control)
    ξ0, U, N = scaling_ξ0(c), c.U, c.N # Definition of the constants
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    return -N * sqrt(U) * b(t) / # Numerator 
           sqrt(
        ξ0^2 * U -          # First term in the denominator
        2 * im * b(t) * db(t)  # Second term in the denominator
    )
end
"""
    `hermite_arg(t, c::Control)`
This function maps the argument of the Hermite polynomial of the STA wave function to the argument of the standard Gaussian function used to evaluate the Fourier transform of the product between a Gaussian function and a Hermite polynomial.
"""
function hermite_arg(t, c::Control)
    ξ0, N = scaling_ξ0(c), c.N # Definition of the constants
    b(t) = auxiliary(t, c) # Auxiliary function
    return ξ0 / (N * b(t))
end

#=
## 2. Exclusively time dependent part

By what I have already described in a different file, we can simplify the product of the time dependent parts of the STA functions.
This simplification can be described as a product of three factors:

1. the normalisation factor, which is the same as the one for the Harmonic Oscillator 
$$\frac{\xi_0}{\sqrt{\pi } \sqrt{b(t)} \sqrt{2^n n! b(t)}}$$

2. the imaginary phase, given by 
$$\exp^\left\{i n \int_0^t \frac{\xi_0^2 U}{2 b(\tau )^2} \, d\tau \right\}$$

3. The factor coming from the Fourier transform 
$$ (-i)^n|\alpha|^2 \left( \alpha^{*2}\beta^{2} -1\right)^{n/2}$$

The first and the third one are easy to implement, the second one is a little bit harder as there is an implicit integral in the exponential.
To overcome this problem, I need to figure out what is the best way to implement the integral.
I think my best take is to try to interpolate the integral function and then define a new one, instead of calling the `quadgk` function recursively.

=#

phase_integrand(t, c::Control) = 1 / auxiliary(t, c)^2
using QuadGK
int(t) = quadgk(t -> phase_integrand(t, c), 0, t)[1]
int(0.2)

