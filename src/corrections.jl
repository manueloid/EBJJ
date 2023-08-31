#=
# Calculation of the corrections

In the other files I have defined all the functions that are going into the calculations of the corrections.
Here I will only define the integrals and other functions that I am going to use to calculate the eSTA corrections.

I am going to follow the same steps I followed in the papers I wrote.

In particular, I need to basically calculate two types of integrals:

1. $ G_n = \int_0^{t_f} dt \langle \chi_n | \Delta H | \chi_0 \rangle $
2. $ K_n = \int_0^{t_f} dt \langle \chi_n | \Nabla H | \chi_0 \rangle $

The definitions of the integrals will be given in the relevant sections.
=#
#=
### 1. $ G_n $

For this kind of integral, we need to calculate the difference between the two Hamiltonians of the system, the original Hamiltonian and the approximated one.

$$
    \Delta H = H_{2} - H_{0} =  - 2 J(t)\left[ e^{-i\hat{p}} b_{h}(z) + b_{h}(z)e^{i\hat{p}} - \partial_{z}^{2}\right] 
    \tag{1}
$$
=#
#=
### 2. $ K_n $

To calculate this term, we need to evaluate the gradient with respect to the change in the control parameters.
Assuming the correction to the control function is a polynomial of order $ \nu $, we can write the $ n $-th term of the gradient as:
$$
	\partial_{\lambda_{j}} H_{2}  =
	-2\prod_{\substack{k=1 \\ k\neq j }}^{\nu}\frac{t-t_{k}}{t_{j} - t_{k}} % Lagrange polynomial term 
	\left[ e^{-i\hat{p}} b_{h}(z) + b_{h}(z)e^{i\hat{p}} \right] % term depending on the b function
    \tag{2}
$$
=#
#=
### 3. Code implementation

By looking at equations (1) and (2), we can see that there are some terms that are common to both integrals.
In particular, the term $ e^{-i\hat{p}} b_{h}(z) + b_{h}(z)e^{i\hat{p}} $ appears in both $ G_n $ and $ K_n $.

Then I need to calculate the part of the gradient of the control function, using the Lagrange polynomial.
Assuming a polynomial of order $ \nu $, the $ i $-th term of the gradient is nothing more than the $ i $-th Lagrange polynomial term over the interval $ [0, t_f] $.

In the following code, I am going to define three functions:

1. `bh_integrand`: this function will return the part of the integrand depending on the $ b_{h}(z) $ function, this term is common to both $ G_n $ and $ K_n $.
2. `sd_integrand`: this function will return the part of the integrand depending on the $ \partial_{z}^{2} $ term, this term is only present in $ G_n $.
3. `gradient`: this function will return the part of the integrand depending on the gradient of the control function, this term is only present in $ K_n $.

The overall plan is to only define the rigth hand side of the integrand, i.e. the one on which the Hamiltonian acts upon.
This is because the left hand side and the time dependent part are in common to the two integrals to calculate the corrections.
I think though that it does not make sense to define general functions.
I would rather define the $ G_n $ and $ K_n $ integrals in place, so that I do not have to define a lot of functions that I am not going to use.
=#

bh(z::Float64, h::Float64) = abs(z) <= 1 ? √((1 + z + h) * (1 - z)) : 0.0 # One-liner function that returns the piecewise function bₕ(z)

"""
    gradient(tarray::Array{Float64, 1})
Calculate the gradient of the control function, given an array of points in the interval  `0, tf`.
The output is an array of anonymous functions, where each function is the  i -th term of the gradient.
This is basically an array of Lagrange polynomials that are 1 at the i-th point and 0 at all the others.
"""
function gradient_int(tarray::Array{Float64,1})
    ncoeffs = length(tarray) - 2
    gradient = Array{Function}(undef, ncoeffs)
    for i in 1:ncoeffs
        yarr = zeros(ncoeffs + 2)
        yarr[i+1] = 1.0
        gradient[i] = t::Float64 -> Lagrange(tarray, yarr)(t)
    end
    return gradient
end

#=
### Definition of the wave function

Here I will define the whole wave function, with no simplifications made and I am going to use it to calculate the corrections.
It will not be the best in terms of code optimization but I want to check if there is any difference in the two approaches.
=#
normalisation(η::ComplexF64, n::Int64) = (real(η) / pi)^(1 / 2) / sqrt(2^n * factorial(n)) * (-im)^n * sqrt(conj(η) / η)^n / abs(η)
gaussian(z::Float64, η::ComplexF64) = exp(-z^2 / (2 * η))
hermite(n::Int64, z::Float64, η::ComplexF64) = SpecialPolynomials.basis(Hermite, n)(z * sqrt(real(η)) / abs(η))
spatial_fourier(n::Int64, z::Float64, η::ComplexF64) = normalisation(η, n) * gaussian(z, η) * hermite(n, z, η)
"""
    simplified(n::Int64, t, x, c::Control)
Return the value of the product between complex conjugate of a STA wave function and the ground state, assuming all the simplifications have been carried out.
"""
function simplified(n::Int64, t, x, c::Control)
    J0, N, U = c.J0, c.N, c.U
    # Definition of the auxiliary functions
    b(t) = auxiliary(t, c)
    db(t) = ForwardDiff.derivative(b, t)
    α(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) - 2im * db(t) / (U * b(t))  # Parameter of the Gaussian term
    αc(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) + 2im * db(t) / (U * b(t)) # Complex Conjugate of the Parameter of the Gaussian term 
    imag_phase_integrand(t::Float64) = sqrt(2J0 * N * U) / b(t)^2
    φ(t::Float64) = quadgk(τ -> imag_phase_integrand(τ), 0.0, t)[1]
    return exp(im * n * φ(t)) * spatial_fourier(n, x, αc(t)) * gaussian(x, αc(t))
end

function corrections(k::Vector{ComplexF64}, g::ComplexF64)
    v = real(conj(g) * k)
    hessian = k * k' |> real
    num = v * norm(v)^2
    den = v' * hessian * v
    return num / den
end
corrections(g::ComplexF64, k::Vector{ComplexF64}) = corrections(k, g)
function corrections(n::Int64, c::Control, λs::Int64=5)
    J0, N, U = c.J0, c.N, c.U
    h = 2.0 / N
    # Definition of the auxiliary functions
    b(t) = auxiliary(t, c)
    db(t) = ForwardDiff.derivative(b, t)
    ddb(t) = ForwardDiff.derivative(db, t)
    J = control_function(b, ddb, c)
    gradient_functions = gradient_int(collect(0.0:c.T/(λs+1):c.T))
    grad(t::Float64) = [g(t) for g in gradient_functions]
    α(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) - 2im * db(t) / (U * b(t))  # Parameter of the Gaussian term
    αc(t::Float64) = 2 / b(t)^2 * sqrt(2J0 * N / U) + 2im * db(t) / (U * b(t)) # Complex Conjugate of the Parameter of the Gaussian term 
    imag_phase_integrand(t::Float64) = sqrt(2J0 * N * U) / b(t)^2
    φ(t::Float64) = quadgk(τ -> imag_phase_integrand(τ), 0.0, t)[1]
    lhs(t, z) = exp(im * n * φ(t)) * spatial_fourier(n, z, αc(t)) * gaussian(z, αc(t)) # left hand side of the integrand, same for both
    rhs_g(t, z) = -2 * J(t) * (
                      bh(z, h) * gaussian(z + h, α(t)) + bh(z - h, h) * gaussian(z - h, α(t)) -
                      h^2 * sd_groundstate(z, α(t))
                  )
    rhs_k(t, z) = -2 * grad(t) * (bh(z, h) * gaussian(z + h, α(t)) + bh(z - h, h) * gaussian(z - h, α(t)))
    gn = hcubature(var -> lhs(var[1], var[2]) * rhs_g(var[1], var[2]), [0.0, -1.0e1], [c.T, 1.0e1], atol=1e-7)[1]
    kn = hcubature(var -> lhs(var[1], var[2]) * rhs_k(var[1], var[2]), [0.0, -1.0e1], [c.T, 1.0e1], atol=1e-7)[1]
    return corrections(gn, kn)
end
function corrections(narr::Vector{Int64}, c::Control; λs::Int64=5)
    corrs = zeros(λs)
    for n in narr
        corrs += corrections(n, c)
    end
    return corrs
end


#=
#### 4.2. $ K_n $

The code to calculate the $ K_n $ is very similar to the one used to calculate the $ G_n $, the only big difference is that I have to repeat the calculation for the different terms of the gradient.

The nice part is that I do not have to use the function depending on the second derivative.

Similarly to what I have done with the $ G_n $ term, I will define a function that takes a lot of arguments.
=#

#=
### 5. Correction calculations

We have now calculated the integrals $ G_n $ and $ K_n $ and we can now calculate the corrections.

In order to calculate the corrections we need to evaluate two main quantities related to the $ G_n $ and $ K_n $ integrals.
These two quantities are

1. The Hessian matrix, the elements of which are defined as
$$
\mathrm{H}_{l,k} = - \sum_{n = 1}^{\mathcal{N}}\left(  \vec{K}^{*}_{n}\right)_{k}\left(  \vec{K}_{n} \right)_{l}\right
$$

(this formulation is quite different from the one from the paper as we can avoid to consider the second derivative term )

2. The vector term
$$
	\vec{v} = \sum_{n=1}^{\mathcal{N}}\text{Re}(G_{n}^{*}\vec{K}_{n})
$$

In both cases, the value $ \mathcal{N} $ is the number of energy levels we are considering in the calculation of the corrections.

#### 5.1. Code implementation

I think it would make more sense to define all the relevant parameters in one big function, and only then evaluate the values for $ G_n $ and $ K_n $, to then finally evaluate the corrections.
=#

