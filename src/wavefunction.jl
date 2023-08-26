#=
# STA Wave function 

Now I am ready to start all the formulation with only one parameter, which is the complex number $ \eta $.
The transform has been working ok, now I need to define the functions that will be used to evaluate the integrand.

In the [test file](./test/wavefunction_tests/fourier_checks.jl) I have set up all the full wave functions, in case they are needed.
In this file though, I will need only some of them, and I will define them as I go.

I will start from what I called the _lefthand side_ of the integrand.
By taking into account the fact that the integrand we need to calculate is of the form $ \langle \chi_n | \hat{O} | \chi_0 \rangle $, we call _lefthand side_ the function $ \chi_n^{*} $ times the part of $ \chi_0 $ upon which the operator doesn't act. 

For the moment I am not going to introduce the imaginary phase as I do not have clear in mind what is the best way to implement it.

I will also start by defining functions that only take general complex values as input, to make the code a little more flexible.
The whole _lefthand side_ function is given by three major terms:
1. the normalisation term
2. the Gaussian term
3. the Hermite polynomial term

The way these factors are derived is described in other files.
When I will need to evaluate the integrand, I will pass the explicit value of the parameter $ \eta $ and it will be the explicit complex conjugate.
By doing that I will be able to avoid the complex conjugation of the whole function, which is a little more expensive.
=#

normalisation(η::ComplexF64, n::Int64) = (real(η) / pi)^(1 / 2) / sqrt(2^n * factorial(n)) * (-im)^n * sqrt(conj(η) / η)^n / abs(η)
gaussian(z::Float64, η::ComplexF64) = exp(-z^2 / (2 * η))
hermite(n::Int64, z::Float64, η::ComplexF64) = SpecialPolynomials.basis(Hermite, n)(z * sqrt(real(η)) / abs(η))
spatial_fourier(n::Int64, z::Float64, η::ComplexF64) = normalisation(η, n) * gaussian(z, η) * hermite(n, z, η)

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
sd_groundstate(z::Float64, η::ComplexF64) = (z^2 - η) * exp(-z^2 / (2 * η)) / η^2

#=
### 2 Imaginary phase

As said earlier, this function is quite tricky from a numerical point of view, as the implicit integral is not easy to evaluate.
To overcome this problem, I need to figure out what is the best way to implement the integral.
I think my best take is to try to interpolate the integral function and then define a new one, instead of calling the `quadgk` function recursively.
After some tests, I found that the interpolation with 1000 points is the best compromise between speed and accuracy.

Regardless, I think the best thing to do is to define this function when I actually need it, i.e. in the body of the functions in which it is used.

### 3 Example of implementation

In the following I will only set up some examples of implementation of the functions defined above, hopefully it will help in the future.

function test_ground_state(n::Int64, t, c::Control)
    J0, N, U = c.J0, c.N, c.U
    b(t) = auxiliary(t, c)                                                    # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t)                                      # Derivative of the auxiliary function
    α(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) - im * db(t) / (U * b(t))  # Parameter of the Gaussian term
    αc(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) + im * db(t) / (U * b(t)) # Complex Conjugate of the Parameter of the Gaussian term 
    lhs(x) = lhs_spatial_fourier(n, x, αc(t))                                 # Lefthand side of the integrand
    rhs(x) = gaussian(x, α(t))                                                # Righthand side of the integrand
    return quadgk(x -> lhs(x) * rhs(x), -Inf, Inf,  atol=1e-7)[1] |> real     # Integral
end
=#
