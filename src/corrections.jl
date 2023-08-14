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
3. `gradient_integrand`: this function will return the part of the integrand depending on the gradient of the control function, this term is only present in $ K_n $.

Moreover, in order to calculate the actual integral, I need to import use the functions defined in the `wavefunctions.jl` file, in particular I need the `time_dependent` function, which is goiong to be common to both $ G_n $ and $ K_n $.
=#

scaling_ξ0(c::Control) = sqrt(8c.J0 * c.N / c.U)
bh(z::Float64, h::Float64) = abs(z) <= 1 ? √((1 + z + h) * (1 - z)) : 0.0 # One-liner function that returns the piecewise function bₕ(z)
"""
    bh_integrand(n::Int64, z, η::ComplexF64)
Return the value of the integrand depending on the `bh` term evaluated at position `z`, given a shift `h`.
It needs a complex number η as input as it will be easier to pass that value later on as a function of time
"""
function bh_integrand(n::Int64, z::Float64, h::Float64, η::ComplexF64)
    return conj(spatial_fourier(n, z, η)) *   # Left part of the integrand
           (
        bh(z - h, h) * ground_state(z - h, η) +   #Second part of the right term of the integrand
        bh(z, h) * ground_state(z + h, η)  #Second part of the right term of the integrand
    )
end
"""
    gradient(tarray::Array{Float64, 1})
Calculate the gradient of the control function, given an array of points in the interval  `[0, tf]`.
The output is an array of anonymous functions, where each function is the  i -th term of the gradient.
"""
function gradient(tarray::Array{Float64,1})
    ncoeffs = length(tarray) - 2
    gradient = Array{Function}(undef, ncoeffs)
    for i in 1:ncoeffs
        yarr = zeros(ncoeffs + 2)
        yarr[i+1] = 1.0
        lagr(t) = Lagrange(tarray, yarr)(t)
        gradient[i] = t::Float64 -> lagr(t)
    end
    return gradient
end

#=
### 4. Code to calculate the $ G_n $ and $ K_n $ 

In this section I am going to define the functions that will calculate the integrals $ G_n $ and $ K_n $.
This is just a composition of the functions I defined here and in the other files.

In this case I am going to define the functions in such a way that it will take the control parameter `ControlFull`as one of the input, and then pass the corresponding values in the functions.
=#
"""
    G_factor(n::Int64, c::Control; λs::Int64=5)
Return the value of the integral `Gₙ` given the energy level `n` and the control parameter `c`.
"""

