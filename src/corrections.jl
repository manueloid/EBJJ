#=
# Calculation of the corrections

In the other files I have defined all the functions that are going into the calculations of the corrections.
Here I will only define the integrals and other functions that I am going to use to calculate the eSTA corrections.

I am going to follow the same steps I followed in the papers I wrote.

In particular, I need to basically calculate two types of integrals:

1. $ G_n = \int_0^{t_f} dt \langle \chi_n | \Delta H | \chi_0 \rangle $
2. $ K_n = \int_0^{t_f} dt \langle \chi_n | \Nabla H | \chi_0 \rangle $

The definitions of the integrals will be given in the relevant sections.

### 1. $ G_n $

For this kind of integral, we need to calculate the difference between the two Hamiltonians of the system, the original Hamiltonian and the approximated one.

$$
    \Delta H = H_{2} - H_{0} =  - 2 J(t)\left[ e^{-i\hat{p}} b_{h}(z) + b_{h}(z)e^{i\hat{p}} - \partial_{z}^{2}\right] 
    \tag{1}
$$


### 2. $ K_n $

To calculate this term, we need to evaluate the gradient with respect to the change in the control parameters.
Assuming the correction to the control function is a polynomial of order $ \nu $, we can write the $ n $-th term of the gradient as:
$$
	\partial_{\lambda_{j}} H_{2}  =
	-2\prod_{\substack{k=1 \\ k\neq j }}^{\nu}\frac{t-t_{k}}{t_{j} - t_{k}} % Lagrange polynomial term 
	\left[ e^{-i\hat{p}} b_{h}(z) + b_{h}(z)e^{i\hat{p}} \right] % term depending on the b function
    \tag{2}
$$

### 3. Considerations about the code

By looking at equations (1) and (2), we can see that there are some terms that are common to both integrals.
In particular, the term $ e^{-i\hat{p}} b_{h}(z) + b_{h}(z)e^{i\hat{p}} $ appears in both $ G_n $ and $ K_n $.

This is the reason why I implemented the function that follows
=#

bh(z::Float64, h::Float64) = abs(z) <= 1 ? √((1 + z + h) * (1 - z)) : 0.0 # One-liner function that returns the piecewise function bₕ(z)
"""
    `bh_integrand(n::Int64, z, η::ComplexF64)`
Return the value of the integrand depending on the `bh` term evaluated at position `z`, given a shift `h`.
It needs a complex number η as input as it will be easier to pass that value later on as a function of time
"""
function bh_integrand(n::Int64, z::Float64, h::Float64, η::ComplexF64)
    return conj(EBJJ.spatial_fourier(n, z, η)) *   # Left part of the integrand
           (
        bh(z - h, h) * EBJJ.ground_state(z - h, η) +   #Second part of the right term of the integrand
        bh(z, h) * EBJJ.ground_state(z + h, η)  #Second part of the right term of the integrand
    )
end

