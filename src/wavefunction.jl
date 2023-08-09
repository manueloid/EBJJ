#=
# STA Wave function 

I need to things slowly and carefully, and write tests along the way.
In this file I will only focus on defining the wave function, in a different one I will actually evaluate the eSTA corrections.

## 1. Constant factor and helpful functions

First I will start by defining the constant factor $\xi_0$ and the functions $\alpha(t)$ and $\beta(t)$ that go into the Fourier transform of the STA wave function.
The names of the function will be (respectively) `scaling_ξ0`, `gaussian_arg` and `hermite_arg`.
I have the feeling that it will make more sense to define everything in one big function instead of defining the constants and the functions separately, as I think it will save some precompilation time.
It will also make sense as I do not have to define the auxiliary function `b(t)` twice.
I have already found out that in this case there is not much to simplify, so I will just write the integrand as it is.

I am going to use the fact that $ \alpha ^2 $ is the inverse of the argument of the Gaussian part in the STA wave function, and if we are going to call the argument $ \eta^2 $, the parameter $ \beta $ is nothing more than $ \sqrt{\Re{\eta^2}} $.

Thus, we can simplify the term $ \alpha ^2 \beta ^2 -1 $ that appears twice in the Fourier transform of the product of the Gaussian term and the Hermite polynomial as $ \eta^{2*} / \eta^2 $
=#

"""
    `spatial_fourier(n::Int64, t, z, η::ComplexF64)`
Evaluate the analytic solution of the Fourier transform of the STA wave function for a complex number η, for a time `t` and at position `z`.
The value η is nothing more than 1/α² in the notes, while β is the real part of η.
"""
function spatial_fourier(n::Int64, z, η::ComplexF64)
    γ = sqrt(conj(η) / η) # this is the √α²β² - 1 factor
    num = sqrt(real(η)) / η # This is the numerator inside the Hermite polynomial
    return exp(-z^2 / (2 * η)) * SpecialPolynomials.basis(Hermite, n)(z * num / γ)
end
ground_state(z, η::ComplexF64) = spatial_fourier(0, z, η)

#=
### 1.1 Integrand depdending on the $ b_h(z) $ function.

I will start by defining a function that given (among other things) the time `t` and a position `x` returns the value of the integrand 
$$
\psi_{n}^{*}(t,x)   \left[ e^{-i\hat{p}} b_{h}(z) + b_{h}(z)e^{i\hat{p}} \right] \psi_{0}(t,x)
$$
where $ \psi_n(t,z) is the spatial part of the STA wave function.
=#

scaling_ξ0(c::Control) = sqrt(8c.J0 * c.N / c.U)

#=
### 1.2 Integrand depending on the second derivative of the ground state

In this section we are going to evaluate the integral of the form 
$$ \langle \psi_n | \frac{d^2}{dz^2} | \psi_0 \rangle $$.

This integral is non zero if $ n $ is even.

Moreover, I could analytically take the second derivative of the ground state, to get
$$\alpha ^2 e^{-\frac{1}{2} \alpha ^2 x^2}\left( \alpha^2 x^2 - 1\right)$$

we can thus obtain a single function to evaluate the integrand, which is given by
    $$
        \mathcal{H}_{n}^{*}\left(\frac{\alpha^{2}\betax}{\sqrt{2\alpha^{2}\beta^{2} - 1}}\right)
       \alpha^2  \exp{- \frac{x^2}{2} \Re{\alpha^2}} \left( \alpha^2 x^2 - 1\right)
    $$

#### 1.2.1 Code implementation

I am going to write a function that returns the value of this integrand given the time `t` and the position `z`.
I do not think there is any need to split the integral into two addends.

On the other hand though, I am not going to reuse any of the functions I used to calculate the part depending on the $ b_h $ factor, as there are more simplifications I can carry out in this section.

I do not know if there is any need to implement the multiple dispatched versions of the functions, but I am going to leave them there just in case.

# Multiple dispatched version of the previous functions

"""
    `bh_integrand(n::Int64, t, z, c::Control)`
Evaluate the integrand of the Fourier transform of the STA wave function for a time `t` and at position `z` given the energy level `n`
This is the multiple dispatched version of the previous one.
"""
function bh_integrand(n::Int64, t, z, c::Control)
    ξ0, U, h = scaling_ξ0(c), c.U, 1 / c.N # Definition of the constants
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    η(t) = ξ0^2 / b(t)^2 - 2im * db(t) / (U * b(t))
    return bh_integrand(n, z, h, η(t))
end
"""
    `sd_integrand(n::Int64, t, z, c::Control)`
Evaluate the integrand of the Fourier transform of the STA wave function for a time `t` and at position `z` given the energy level `n`, for a control parameter `c`.
This is the multiple dispatched version of the previous one.
"""
function sd_integrand(n::Int64, t, z, c::Control)
    ξ0, U = scaling_ξ0(c), c.U # Definition of the constants
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    η(t) = ξ0^2 / b(t)^2 - 2im * db(t) / (U * b(t))
    return sd_integrand(n, z, η(t))
end
=#

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

I would like to try to implement the code in such a way that no allocations are made, for example I would like to not to define the auxiliary function.
I could do that by either passing a value or a function as an argument.
The second option is a little bit more trickier because at the moment I do not know how to pass a function.
It seems like I only need to use the `::Function` type annotation, so I will keep on going with it.
=#

normalisation(n::Int64, t, ξ0::Float64, b) = ξ0 * (pi * b(t)^2 * 2^n * factorial(n))^-0.5

#=
### 2.2 Imaginary phase
As said earlier, this function is quite tricky from a numerical point of view, as the implicit integral is not easy to evaluate.
To overcome this problem, I need to figure out what is the best way to implement the integral.
I think my best take is to try to interpolate the integral function and then define a new one, instead of calling the `quadgk` function recursively.
After some tests, I found that the interpolation with 1000 points is the best compromise between speed and accuracy.

I will first define the general algorithm to return the integral function given the function to integrate and the time interval.
I decided to have the function returning the interpolated function so that I do not have to interpolate it every time I need to evaluate it.
=#

"""
    `interpolation_integral(tf::Float64, b::Function64; npoints=1000)`
Return a function which is the interpolation of the function f(t) = ∫₀ᵗ b(τ) dτ in the interval [0, tf].
The function takes the argument `npoints` which is the number of points used for the interpolation. The default value is 1000 as it is the best trade-off between speed and accuracy.
"""
function interpolation_integral(tf::Float64, b; npoints=1000)
    trange = range(0.0, tf, length=npoints)
    integral_values = [quadgk(t -> b(t), 0.0, t)[1] for t in trange]
    itp = linear_interpolation(trange, integral_values)
    return t::Float64 -> itp(t)
end
"""
    `imaginary_phase(n, t, c::Control; npoints=1000)`
This function returns the imaginary phase of the product of the nth STA wave function and the ground state, as defined in the notes.

The argument it takes are 
- `n` which is the level of excitation of the wave function
- `t` which is the time at which the wave function is evaluated
- `c` which is the control object, containing all the parameters of the system
- `φ` which is the integral function that goes into the exponential, it could be either the interpolated function or the actual integral function.
"""
function imaginary_phase(n, t, c::Control, φ)
    ξ0, U = scaling_ξ0(c), c.U # Definition of the constants
    return exp(im * n * ξ0^2 * U / 2 * φ(t))
end

#=
### 2.3 Fourier transform factor
This term comes from the Fourier transform of the product between a Gaussian function and a Hermite polynomial.
It is given by (using the $\alpha$ and $\beta$ notation from the notes)
$$ (-i)^n|\alpha|^2 \left( \alpha^{*2}\beta^{2} -1\right)^{n/2}$$.

But if we use the $ \eta $ notation, we can write it as 
    $$ (-i)^n \Re{\eta^2} \left(\frac{\eta^{2}}{\eta^{2*}})^{n/2}$$.

Where we would like to point out that we swapped the term in the fraction as we are considering the complex conjugate of the $\eta$ function.

The function I am going to define will take a general complex number `η` as an argument, as it is easier to implement it this way.
=#
"""
    fourier_factor(n, η::Complex)
Return the Fourier transform factor as a function of the energy level `n` and a general complex number `η`.
"""
fourier_factor(n::Int64, η::Complex) = (-im)^n * real(η) * (η / conj(η))^(n / 2)

#=
## 2.4 Time dependent part 
Here I will just implement the whole time dependent part of the wave function, just to make it a little bit easier to read.

This function will be nothing more than a combination of the previous ones, in which I am going to pass the specific functions and the complex numbers.
A possible implementation will be something like this:
=#

"""
    time_dependent(n::Int64, t, c::Control)
Return the time dependent part of the product between the nth STA wave function and the ground state, as defined in the notes.
"""
function time_dependent(n::Int64, t, c::Control)
    ξ0, U = scaling_ξ0(c), c.U # Definition of the constants
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    η(t) = ξ0^2 / b(t)^2 - 2im * db(t) / (U * b(t))
    φ = interpolation_integral(t, b)
    return normalisation(n, t, ξ0, b) * imaginary_phase(n, t, c, φ) * fourier_factor(n, η(t))
end
