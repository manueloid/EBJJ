#=
# STA Wave function 

Now I am ready to start all the formulation with only one parameter, which is the complex number $ \alpha^2 $.
The transform has been working ok, now I need to define the functions that will be used to evaluate the integrand.

In the [test file](./test/wavefunction_tests/fourier_checks.jl) I have set up all the full wave functions, in case they are needed.
In this file though, I will need only some of them, and I will define them as I go.

The whole STA wave function, result of the Fourier transform of the solution of the Harmonic oscillator in momentum representation, is given by:
$$
    \chi_n(z,t) =
    \left[\frac{R_\alpha^2}{\pi}\right]^{1/4}
    \frac{1}{\sqrt{2^n n!}} % Hermite polynomial normalisation factor   
	\exp\left\{-i(n + 1/2)\ \int_{0}^{t}d \tau ~ U R_\alpha^2(\tau)\right\} % time dependent phase factor
    \frac{i^n}{\sqrt{\hbar}} \frac{1}{\alpha} \left(\frac{\alpha^{2*}}{\alpha^2}\right)^{n/2} % normalisation factor
    \exp\left\{-\frac{z^2}{2 \hbar^2 \alpha^2} \right\} % Gaussian term
    \mathcal{H}_n\left(\frac{R_\alpha}{|\alpha^2|}~\frac{z}{\hbar}\right) 
    \tag{1}
$$
where $ \alpha $ is a complex number and $ R_\alpha $ is its real part.

I will start from what I called the _lefthand side_ of the integrand.
By taking into account the fact that the integrand we need to calculate is of the form $ \langle \chi_n | \hat{O} | \chi_0 \rangle $, we call _lefthand side_ the function $ \chi_n^{*} $ times the part of $ \chi_0 $ upon which the operator doesn't act. 

I will start by defining functions that only take general complex values as input, to make the code a little more flexible.
The whole _lefthand side_ function is given by three major terms:
1. the normalisation terms
2. the Gaussian term
3. the Hermite polynomial term

In the following section I will give an overview of the functions I need.
We will start from the normalisation terms which is the part that is not depending on the position variable $ z $.
=#
#=
## 1 Time dependent term

We know that the operator $ \hat{O} $ acts only on the position variable, then we can simplify the expression by considering only the time dependent part of the wave function and multiplying them together.
In particular, we need to multiply the time dependent part of the wave function with the complex conjugate of the time dependent part of the wave function.

So we can pair them together and obtain one function depending on the time variable $ t $:
$$
    \left[\frac{R_\alpha^2}{\pi}\right]^{1/2}
    \frac{1}{\sqrt{2^n n!}} % Hermite polynomial normalisation factor
    \frac{(-i)^n}{\hbar} \frac{1}{|\alpha|^2}
    \left(\frac{\alpha^{2}}{\alpha^{2*}}\right)^{n/2} % normalisation factor
    \exp\left\{ i n \int_{0}^{t}d \tau ~ U R_\alpha^2(\tau)\right\} % time dependent phase factor
    := \Phi(t)
    \tag{2}
$$

To give a little bit more perspective on the problem, we need to calculate the following integral 
$$
 \langle \chi_n(z,t) | \hat{O} | \chi_0(z,t) \rangle
 \tag{3}
$$
where $ \hat{O} $ is a generic operator acting on the position variable $ z $.

We have already written down the STA wave function $ \chi_n(z,t) $ in eq. (1), now we need to write down both the ground state $ \chi_0(z,t) $ and the complex conjugate $ \langle \chi_n(z,t) | $.

The first one is given by
$$
    \chi_0(z,t) =
    \left[\frac{R_\alpha^2}{\pi}\right]^{1/4}
    \exp\left\{-\frac{i}{2}\ \int_{0}^{t}d\tau ~ U R_\alpha^2(\tau)\right\} % time dependent phase factor
    \frac{1}{\sqrt{\hbar}} \frac{1}{\alpha} % normalisation factor
    \exp\left\{-\frac{z^2}{2 \hbar^2 \alpha^2} \right\} % Gaussian term
    \tag{4}
$$
while the second one is given by

$$
    \chi_n(z,t)^* =
    \left[\frac{R_\alpha^2}{\pi}\right]^{1/4}
    \frac{1}{\sqrt{2^n n!}} % Hermite polynomial normalisation factor   
	\exp\left\{i(n + 1/2)\ \int_{0}^{t}d \tau   ~ U R_\alpha^2(\tau)\right\} % time dependent phase factor
    \frac{(-i)^n}{\sqrt{\hbar}} \frac{1}{\alpha^*} \left(\frac{\alpha^{2}}{\alpha^{2*}}\right)^{n/2} % normalisation factor
    \exp\left\{-\frac{z^2}{2 \hbar^2 \alpha^{2*}} \right\} % Gaussian term
    \mathcal{H}_n\left(\frac{R_\alpha}{|\alpha^2|}~\frac{z}{\hbar}\right) 
    \tag{5}
$$
Where $ \hbar $ in this case is $ 1/N$.

Finally, the whole time dependent term is given by
$$
    \left[\frac{R_\alpha^2}{\pi}\right]^{1/2}
    \frac{1}{\sqrt{2^n n!}} % Hermite polynomial normalisation factor
    \frac{(-i)^n}{\hbar} \frac{1}{|\alpha|^2}
    \left(\frac{\alpha^{2}}{\alpha^{2*}}\right)^{n/2} % normalisation factor
    \exp\left\{ i n \int_{0}^{t}d \tau ~ U R_\alpha^2(\tau)\right\} % time dependent phase factor
    := \Phi(t)
    \tag{6}
$$

The types of integrals we need to calculate are then of the form $\langle \chi_n(z,t) | \hat{O} | \chi_0(z,t) \rangle$.
They are given by
$$
    \int_{0}^{t_f} ~ dt \Phi(t) \langle \chi_n(z,t) | \hat{O} | \chi_0(z,t) \rangle = 
    \int_{0}^{t_f} ~ dt \Phi(t) \int_{-\infty}^{\infty} dz ~ \mathcal{H}_n\left(\frac{R_\alpha}{|\alpha^2|} \frac{z}{\hbar}\right) 
    \exp\left\{-\frac{z^2}{2 \hbar^2 \alpha^{2*}} \right\} % Gaussian term conjugate
    \hat{O} \exp\left\{-\frac{z^2}{2 \hbar^2 \alpha^2} \right\} % Gaussian term
    \tag{7}
$$
=#
#=
### 1.1 Imaginary phase factor

I am going to define this part of the wave function inside the function that evaluates the integrand.
I will define a function that returns an integral function that will go inside the integral later on.

### 1.2 Wrap up

To conclude, for the left-hand side, we need to define the time dependent term of the product between the ground state and the general complex conjugate of $ \chi_n $.
For the Hermitian polynomial part of $ \langle \chi_n | $, we do not need any complex conjugate as the argument is purely real.
For the Gaussian part of $ \langle \chi_n | $, we need to define the complex conjugate of the Gaussian term.
I could either take the complex conjugate of the whole term or directly pass the complex conjugate of the parameter $ \alpha^2 $.

The second option is the more viable as it less computationally expensive.
=#
#=
### 1.3 Code implementation
I will start by defining a function that takes a complex number `α2` and its complex conjugate `α2c` as arguments and returns the time dependent term of the wave function apart from the imaginary phase factor, as I will define that later.

=#

he(n::Int64, z::Float64) = SpecialPolynomials.basis(Hermite, n)(z) # one liner to speed up the writing of the polynomial
"""
    norm(n::Int64, α2::ComplexF64, α2c::ComplexF64, h::Float64)
Return the time dependent part of the product ⟨χₙ|χ₀⟩, for a general complex number `α2` its complex conjugate `α2c` and the energy level `n`.
The function also takes the parameter `h`, for later use.
The function takes two complex numbers as input to reduce the computational cost.
"""
function norm(n::Int64, α2::ComplexF64, α2c::ComplexF64, h::Float64)
    r = sqrt(real(α2))
    return (r^2 / pi)^(1 / 2) * (1 / sqrt(2^n * factorial(n))) * ((-im)^n / h) * (1 / abs(α2)) * (α2 / α2c)^(n / 2)
end
"""
    ga_wh(z::Float64, α2::ComplexF64)
Return the Gaussian part of the wave function at a given position for a general complex number `α2`.
The idea is to pass the variable `z/h` when calling the function, so there is no need to pass the value `h` in this functio
For the complex conjugate of this function, is just enough to pass the complex conjugate of the parameter `α2`.
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

#=
## 2 Integrand depending on the second derivative of the ground state

In this section we are going to evaluate the integral of the form 
$$ \langle \psi_n | \frac{d^2}{dz^2} | \psi_0 \rangle $$.

This integral is non zero if $ n $ is even.

Moreover, I could analytically take the second derivative of the ground state.
The only take we need to be careful of, is the fact that in this case the exponential is of the form $$ \exp\left\{-\frac{z^2}{2 \hbar^2 \alpha^{2}} \right\} $$, so there is an extra $\hbar $ term we need to take account of.((z - h*α)*(z + h*α))/(E^(z^2/(2*h^2*α^2))*h^4*α^4)
The resulting function is given by
$$
    \partial_z^2 
    \exp\left\{-\frac{z^2}{2 \hbar^2 \alpha^2} \right\} % Gaussian term
    = \frac{e^{-\frac{z^2}{2 \alpha ^2 h^2}} \left(z^2-\alpha ^2 h^2\right)}{\alpha ^4 h^4}
\tag{8}
$$

#### 1.2.1 Code implementation

Here I am going to only define the function that returns the value of the second derivative of the ground state, I will not define the whole integrand function as I will do that in the [relevant file](./src/corrections.jl).
The function will take only a complex number `η` and a position `z` as arguments, and it will return a complex number.

=#
sd_groundstate(z::Float64, α2::ComplexF64, h::Float64) = (z^2 - α2 * h^2) * exp(-z^2 / (2 * α2 * h^2)) / (α2^2 * h^4)
