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
$$ \xi_0 \left( \pi b^2(t) 2^n n! \right)^{-1/2}$$

2. the imaginary phase, given by 
$$\exp^\left\{i n \int_0^t \frac{\xi_0^2 U}{2 b(\tau )^2} \, d\tau \right\}$$

3. The factor coming from the Fourier transform 
$$ (-i)^n|\alpha|^2 \left( \alpha^{*2}\beta^{2} -1\right)^{n/2}$$

The first and the third one are easy to implement, the second one is a little bit harder as there is an implicit integral in the exponential.

### 2.1 Normalisation factor
The normalisation factor depends on the energy level $ n $ and it will be a function of time.
=#
"""
    `normalisation(n, t, c::Control)` 
Returns the normalisation factor as a function of the energy level `n` and at the time `t`.
"""
function normalisation(n, t, c::Control)
    ξ0, N = scaling_ξ0(c), c.N # Definition of the constants
    b(t) = auxiliary(t, c) # Auxiliary function
    return ξ0 * (pi * b(t)^2 * 2^n * factorial(n))^-0.5
end

#=
### 2.2 Imaginary phase
As said earlier, this function is quite tricky from a numerical point of view, as the implicit integral is not easy to evaluate.
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
It returns a function that can be evaluated at any point in the interval [0, T].
This is only the integral of 1/b(t)^2, so it is not the phase function itself.
"""
function interpolation_integral(c::Control; npoints=1000)
    phase_integrand(t, c::Control) = 1 / auxiliary(t, c)^2
    int(t) = quadgk(t -> phase_integrand(t, c), 0, t)[1]
    trange = range(0.0, c.T, length=npoints)
    integral_values = [int(t) for t in trange]
    itp = linear_interpolation(trange, integral_values)
    return t -> itp(t)
end
"""
    `imaginary_phase_pr(n, t, c::Control; npoints=1000)`
This function returns the imaginary phase of the product of the nth STA wave function and the ground state, as defined in the notes.

The argument it takes are 
- `n` which is the level of excitation of the wave function
- `t` which is the time at which the wave function is evaluated
- `c` which is the control object, containing all the parameters of the system
- `integral_func` which is the integral function that goes into the exponential, it could be either the interpolated function or the actual integral function.

"""
function imaginary_phase_pr(n, t, c::Control, integral_func::Function)
    ξ0, U = scaling_ξ0(c), c.U # Definition of the constants
    return exp(im * n * ξ0^2 * U / 2 * integral_func(t))
end

#=
### 2.3 Fourier transform factor
This term comes from the Fourier transform of the product between a Gaussian function and a Hermite polynomial.
It is a function of the time `t` and of the control object `c`, as well as of the energy level `n`.
=#

"""
    `fourier_factor(n, t, c::Control)`
Return the Fourier transform factor as a function of the energy level `n`, the time `t` and the control object `c`.
"""
function fourier_factor(n, t, c::Control)
    αc(t), β(t) = conj(gaussian_arg(t, c)), hermite_arg(t, c) # Arguments of the Gaussian and Hermite functions
    return (-im)^n * abs(αc(t))^2 * (αc(t)^2 * β(t)^2 - 1)^(n / 2)
end

#=
## 2.4 Time dependent part 
Here I will just implement the whole time dependent part of the wave function, just to make it a little bit easier to read.
=#

"""
    `time_dependent(n, t, c::Control)`
Return the time dependent part of the STA wave function, as a function of the energy level `n` and the time `t`, for a given set of control parameter `c`.
The keyword argument `npoints` is the number of points used to interpolate the integral in the imaginary phase function.
"""
function time_dependent(n, t, c::Control, integral_func::Function)
    return normalisation(n, t, c) * imaginary_phase_pr(n, t, c, integral_func) * fourier_factor(n, t, c)
end
precompile(time_dependent, (Int64, Float64, ControlFull, Function))
precompile(time_dependent, (Int64, Float64, ControlInt, Function))

#=
## 3. Spatial part of the wave function
In this case, the spatial part of the wave function is given by the Hermite polynomial and the Gaussian term.
The first one depends on the energy level `n`, time and space `t` and `z`, while the second one only depends on the time `t` and the space `z`.
Let us call this function $f_n(t,z) = g(t, z) h_n(t, z)$, where $g(t, z)$ is the Gaussian term and $h_n(t, z)$ is the Hermite polynomial.
These two functions are defined as:

1. Gaussian term $ g(t,z) =	\exp\left\{-\frac{x^{2}\alpha^{2}}{2}\right\}$
2. Hermite polynomial term $ h_n(t,z) = 
$\mathcal{H}_{n}\left(\frac{\alpha^{2}\betax}{\sqrt{2\alpha^{2}\beta^{2} - 1}}\right)$

What we need to is to implement a function that returns the complex conjugate of $f_n(t,z)$, as this is the term that appears in the integral, and also define a function that only returns $f_0(t,z)$ which - as we already saw - is nothing more than a Gaussian function.
=#
"""
ground_state(t,z,c::Control)
Return the ground state wave function as a function of time `t` and space `z`, for a given set of control parameters `c`.
This function will be then reused to define the complex conjugate of the Gaussian term for the spatial factor of the general STA wave function.
"""
function ground_state(t, z, c::Control)
    α(t) = gaussian_arg(t, c) # Argument of the Gaussian function
    return exp(-z^2α(t)^2 / 2)
end

#=
We can now implement the spatial part of the general STA wave function.
We need to use the `ground_state` function defined above, as well as the `SpecialPolynomials` package to compute the Hermite polynomial.
=#
"""
`hermite(n, t, z, c::Control)`
Return the Fourier antitrasform of the Hermite polynomial of order `n`, as a function of time `t`, space `z` and control parameters `c`.
The function `γ(t)` is the argument of the Fourier transform of the Hermite polynomial.
"""
function hermite(n, t, z, c::Control)
    αc(t), β(t) = conj(gaussian_arg(t, c)), hermite_arg(t, c) # Arguments of the Gaussian and Hermite functions   
    γ(t) = αc(t)^2 * β(t) / sqrt(2 * αc(t)^2 * β(t)^2 - 1) # Argument of the Hermite polynomial  
    return SpecialPolynomials.basis(Hermite, n)(z * γ(t))
end

#=
Now we only need to put the two together to obtain the spatial term of the STA wave function.
=#
"""
`spatial(n, t, z, c::Control)`
Return the spatial part of the STA wave function, as a function of the energy level `n`, time `t`, space `z` and control parameters `c`.
This is already the complex conjugate, as it is the left term of the inner product.
"""
function spatial(n, t, z, c::Control)
    return conj(ground_state(t, z, c)) * hermite(n, t, z, c)
end
precompile(spatial, (Int64, Float64, Float64, ControlFull))
precompile(spatial, (Int64, Float64, Float64, ControlInt))
