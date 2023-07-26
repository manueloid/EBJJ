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
After some tests, I found that the interpolation with 1000 points is the best compromise between speed and accuracy.

I will hence define the function `interpolation_integral(c::Control)` that returns a function that can be evaluated at any point in the interval [0, T].
Then I will define the actual phase function `imaginary_phase(t, c::Control)` that will be used in the wave function.
=#

"""
    `interpolation_integral(c::Control; npoints=1000)`
This function interpolates the values of the integral of the phase function, and returns a function that can be evaluated at any point in the interval [0, T].
This function takes a keyword argument `npoints` which is the number of points used to interpolate the integral, default value is 1000 which is a good compromise between speed and accuracy.
This is only the integral of 1/b(t)^2, so it is not the phase function itself.
"""
function interpolation_integral(c::Control; npoints=1000)
    phase_integrand(t, c::Control) = 1 / auxiliary(t, c)^2
    int(t) = quadgk(t -> phase_integrand(t, c), 0, t)[1]
    trange = range(0.0, c.T, length=npoints)
    integral_values = [int(t) for t in trange]
    itp = linear_interpolation(trange, integral_values)
    return itp
end
"""
    `imaginary_phase(t, c::Control; npoints=1000)`
This function returns the imaginary phase of the STA wave function, as defined in the notes.
The argument it takes are 
- `n` which is the level of excitation of the wave function
- `t` which is the time at which the wave function is evaluated
- `c` which is the control object, containing all the parameters of the system

This function also takes a keyword argument `npoints` which is the number of points used to interpolate the integral.
"""
function imaginary_phase(n, t, c::Control; npoints=1000)
    ξ0, U = scaling_ξ0(c), c.U # Definition of the constants
    itp = interpolation_integral(c, npoints=npoints) # Interpolation of the integral
    return exp(im * n * ξ0^2 * U / 2 * itp(t))
end
