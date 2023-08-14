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

By what I have already described in a [different file](~/Repos/ExternalBJJ/Documents/Formulae/latex_build/formulae.pdf), we can simplify the product of the time dependent parts of the STA functions.
This is

$$
    \sqrt{\frac{\Re[\eta^2]}{\pi}} % normalisation not depending on the excitation
	\frac{-i^n \gamma^{n}}{\sqrt{2^{n}n!}} % normalisation term depending on the excitation
	\exp\left\{i n \int_{0}^{t} \frac{\xi_{0}^{2} U }{2 b(\tau)^{2}} d\tau \right\} % imaginary time evolution factor
    tag{1}
$$ 

where we have
$$
	\eta^2 = \frac{\xi_0^2}{b(t)^2}-\frac{2 i b'(t)}{U b(t)}
    \tag{2}
$$

and $ \gamma $ is $\frac{\eta^{2*}}{\eta^{2}}$.

### 2.1 Code implementation

I found out that in a very high level implementation of the code, I can define both the time dependent part of the product of the wave function by passing just three arguments:
1. the energy level `n`
2. a complex number `η`
3. a real number `k`

The real number `k` is the value of the integral function at the time `t`, and it is used to evaluate the imaginary phase.
To be verbose, the value of `k` is given by (at time `t`)
$$ \int_0^t \frac{\xi_0^2 U}{2 b(\tau )^2} \, d\tau $$

While the value of `η` is given by (at time `t` as well )
$$ \frac{\xi_0^2}{b(t)^2} - \frac{2i}{U} \frac{db(t)}{b(t)} $$

From equation 1, we can see that the normalisation factor can be split into two parts, one where the excitation energy is taken into account, and one where it is not.
I can thus use the multiple dispatch where I can pass only one argument, and one where I can pass both the energy level and the complex numbers.

I need to remember that the name of the variable `\eta` can be a little misleading, as it actually not the complex number $\eta$ but its square.
Moreover, I will define the second normalisation function as taking `γ` as an argument, as it makes the code more readable.

=#

# First part of the normalisation
normalisation(reη::Float64) = sqrt(reη) / sqrt(pi)
normalisation(η::ComplexF64) = normalisation(real(η))

# Second part of the normalisation, the one with the excitation energy
normalisation(n::Int64, γ::ComplexF64) = (-im)^n * (γ)^n / sqrt(2^n * factorial(n))
normalisation(γ::ComplexF64, n::Int64) = normalisation(n, γ)

# Imaginary phase
imaginary_phase(n::Int64, k::Float64) = exp(im * n * k)
imaginary_phase(k::Float64, n::Int64) = imaginary_phase(n, k)

#=
#### 2.1.1 Imaginary phase
As said earlier, this function is quite tricky from a numerical point of view, as the implicit integral is not easy to evaluate.
To overcome this problem, I need to figure out what is the best way to implement the integral.
I think my best take is to try to interpolate the integral function and then define a new one, instead of calling the `quadgk` function recursively.
After some tests, I found that the interpolation with 1000 points is the best compromise between speed and accuracy.

I will first define the general algorithm to return the integral function given the function to integrate and the time interval.
I decided to have the function returning the interpolated function so that I do not have to interpolate it every time I need to evaluate it.

Finally, I chose not to define the final imaginary phase function, I will define it later in the `time_dependent` function.
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
    return itp
end
