#=
# STA Wave function 

I need to things slowly and carefully, and write tests along the way.
In this file I will only focus on defining the wave function and the functions that will go into the integrals, in a different one I will actually evaluate the eSTA corrections i.e. the integrals.

I am going to use the fact that $ \alpha ^2 $ is the inverse of the argument of the Gaussian part in the STA wave function, and if we are going to call the argument $ \eta^2 $, the parameter $ \beta $ is nothing more than $ \sqrt{\Re{\eta^2}} $.

Thus, we can simplify the term $ \alpha ^2 \beta ^2 -1 $ that appears twice in the Fourier transform of the product of the Gaussian term and the Hermite polynomial as $ \eta^{2*} / \eta^2 $

I will first focus on the spatial term of the wave function and only then I will move on to the time dependent part, which is common to all the integrals and can thus be simplified greatly.

Throughout this document, I will not define any function depeding on the actual values of the control parameters, I will only define functions that take general complex and real numbers as arguments.
I will specialise the functions in the [relevant file](./src/corrections.jl), where I will also define the integrals.

## 1. Spatial part

For this term I will only need the product of the Gaussian term and the Hermite polynomial, where the only _tricky_ part is the argument of the Hermite polynomial, which is given by

$$ \frac{\alpha^2 x^2}{\sqrt{2\alpha^2 \beta^2 - 1}} $$

that can be rewritten with the formulation of $ \eta $ as

$$ \frac{\sqrt{\Re[\eta^2]}}{\eta} \frac{1}{\gamma}  $$

where $ \gamma = \sqrt{\eta^{2*} / \eta^2} $.

### 1.1 General spatial part of the wave function and ground state

Here I will define the general spatial part of the ground state of the STA wave function, as well as the ground state for it.
The former is the lefthand side of the integral $ \langle \psi_n |\hat{H} | \psi_0 \rangle $, while the latter is the righthand side and it is only the gaussian term, the Hermite polynomial being one.

#### 1.1.1 Code implementation

I think it makes sense to define a function that takes only the energy level `n` and the position `z` as arguments as well as a general complex number `η`.
I will then define the multiple dispatch versions of the function, where I can pass the arguments in any order.

I would like to point out that `spatial_fourier` is not yet the complex conjugate of the spatial part, and it has to be implemented in due time.
=#

"""
    spatial_fourier(n::Int64, z, η::ComplexF64)
Evaluate the analytic solution of the Fourier transform of the STA wave function for a complex number η, for a time `t` and at position `z`.
The value η is nothing more than 1/α² in the notes, while β is the real part of η.
"""
function spatial_fourier(n::Int64, z::Float64, η::ComplexF64)
    γ = sqrt(conj(η) / η) # this is the √α²β² - 1 factor
    num = sqrt(real(η)) / η # This is the numerator inside the Hermite polynomial
    return exp(-z^2 / (2 * η)) * SpecialPolynomials.basis(Hermite, n)(z * num / γ)
end

# Multiple dispatch version of the previous function
spatial_fourier(n::Int64, η::ComplexF64, z::Float64) = spatial_fourier(n, z, η)
spatial_fourier(η::ComplexF64, n::Int64, z::Float64) = spatial_fourier(n, z, η)
spatial_fourier(η::ComplexF64, z::Float64, n::Int64) = spatial_fourier(n, z, η)
spatial_fourier(z::Float64, η::ComplexF64, n::Int64) = spatial_fourier(n, z, η)
spatial_fourier(z::Float64, n::Int64, η::ComplexF64) = spatial_fourier(n, z, η)

# Ground state of the STA wave function multiply dispatched
ground_state(z::Float64, η::ComplexF64) = exp(-z^2 / (2 * η))
ground_state(η::ComplexF64, z::Float64) = ground_state(z, η)

#=
Since the part of the integrand depending on the $ b_h $ factor does not have any simplification that can be carried out, there is no need to define it here.
This will be done in the [relevant corrections file ](./src/corrections.jl).

### 1.2 Integrand depending on the second derivative of the ground state

In this section we are going to evaluate the integral of the form 
$$ \langle \psi_n | \frac{d^2}{dz^2} | \psi_0 \rangle $$.

This integral is non zero if $ n $ is even.

Moreover, I could analytically take the second derivative of the ground state, to get
$$
    \alpha ^2 e^{-\frac{1}{2} \alpha ^2 x^2}\left( \alpha^2 x^2 - 1\right)
$$

which - in the $ \eta $ formulation - is given by
$$
    \frac{e^{-\frac{z^2}{2 \eta ^2}} \left(z^2-\eta ^2\right)}{\eta ^4}
$$

we can thus obtain a single function to evaluate the integrand, which is given by
$$
    \mathcal{H}_{n}^{*}\left(\frac{\alpha^{2}\beta x}{\sqrt{2\alpha^{2}\beta^{2} - 1}}\right)
    \alpha^2  \exp{- \frac{x^2}{2} \Re{\alpha^2}} \left( \alpha^2 x^2 - 1\right)
$$

#### 1.2.1 Code implementation

Here I am going to only define the function that returns the value of the second derivative of the ground state, I will not define the whole integrand function as I will do that in the [relevant file](./src/corrections.jl).
The function will take only a complex number `η` and a position `z` as arguments, and it will return a complex number.

=#
sd_groundstate(z::Float64, η::ComplexF64) = (z^2 - η^2) * exp(-z^2 / (2 * η^2)) / η^4
sd_groundstate(η::ComplexF64, z::Float64) = sd_groundstate(z, η)

#=
## 2. Exclusively time dependent part

By what I have already described in a [different file](~/Repos/ExternalBJJ/Documents/Formulae/latex_build/formulae.pdf), we can simplify the product of the time dependent parts of the STA functions.
This is

$$
\sqrt{\frac{\Re[\eta^2]}{\pi}} \frac{1}{\Re[\eta^2]} % normalisation not depending on the excitation
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

Again, the remember that 
=#

# First part of the normalisation
normalisation(reη::Float64) = sqrt(reη) / sqrt(pi) / reη^2
normalisation(η::ComplexF64) = normalisation(real(η))

# Second part of the normalisation, the one with the excitation energy
normalisation(n::Int64, γ::ComplexF64) = (-im)^n * (γ)^n / sqrt(2^n * factorial(n))
normalisation(γ::ComplexF64, n::Int64) = normalisation(n, γ)

# Imaginary phase
imaginary_phase(n::Int64, k::Float64) = exp(im * n * k)
imaginary_phase(k::Float64, n::Int64) = imaginary_phase(n, k)

# Full function
function time_dependent(n::Int64, η::ComplexF64, k::Float64)
    γ = sqrt(conj(η) / η)
    return normalisation(η) * normalisation(n, γ) * imaginary_phase(n, k)
end

# Multiple dispatch combinations
time_dependent(n::Int64, k::Float64, η::ComplexF64) = time_dependent(n, η, k)
time_dependent(k::Float64, η::ComplexF64, n::Int64) = time_dependent(n, η, k)
time_dependent(k::Float64, n::Int64, η::ComplexF64) = time_dependent(n, η, k)
time_dependent(η::ComplexF64, k::Float64, n::Int64) = time_dependent(n, η, k)
time_dependent(η::ComplexF64, n::Int64, k::Float64) = time_dependent(n, η, k)

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
    interpolation_integral(tf::Float64, b::Function64; npoints=1000)
Return a function which is the interpolation of the function f(t) = ∫₀ᵗ b(τ) dτ in the interval [0, tf].
The function takes the argument `npoints` which is the number of points used for the interpolation. The default value is 1000 as it is the best trade-off between speed and accuracy.
"""
function interpolation_integral(tf::Float64, b; npoints=1000)
    trange = range(0.0, tf, length=npoints)
    integral_values = [quadgk(t -> b(t), 0.0, t)[1] for t in trange]
    itp = linear_interpolation(trange, integral_values)
    return itp
end
interpolation_integral(c::Control, b; npoints=1000) = interpolation_integral(c.T, b, npoints=npoints)
function interpolation_integral(c::Control; npoints=1000)
    ξ0, U = EBJJ.scaling_ξ0(c), c.U
    b(t::Float64) = auxiliary(t, c)
    to_int(t::Float64) = ξ0^2 * U / b(t)^2
    return interpolation_integral(c, to_int; npoints=npoints)
end

