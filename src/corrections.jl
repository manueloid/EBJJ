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

The overall plan is to only define the rigth hand side of the integrand, i.e. the one on which the Hamiltonian acts upon.
This is because the left hand side and the time dependent part are in common to the two integrals to calculate the corrections.
=#

scaling_ξ0(c::Control) = sqrt(8c.J0 * c.N / c.U)
bh(z::Float64, h::Float64) = abs(z) <= 1 ? √((1 + z + h) * (1 - z)) : 0.0 # One-liner function that returns the piecewise function bₕ(z)
"""
    bh_integrand(n::Int64, z, η::ComplexF64)
Return the value of the integrand depending on the `bh` term evaluated at position `z`, given a shift `h`.
It needs a complex number η as input as it will be easier to pass that value later on as a function of time
"""
function bh_integrand(z::Float64, h::Float64, η::ComplexF64)
    return bh(z - h, h) * ground_state(z - h, η) + bh(z, h) * ground_state(z + h, η)

end
# no need to define the part depending on the second derivative of the ground state as this is already defined elsewhere.

"""
    gradient(tarray::Array{Float64, 1})
Calculate the gradient of the control function, given an array of points in the interval  `0, tf`.
The output is an array of anonymous functions, where each function is the  i -th term of the gradient.
This is basically an array of Lagrange polynomials that are 1 at the i-th point and 0 at all the others.
"""
function gradient(tarray::Array{Float64,1})
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
### 4. Code to calculate the $ G_n $ and $ K_n $ 

In this section I am going to define the functions that will calculate the integrals $ G_n $ and $ K_n $.
This is just a composition of the functions I defined here and in the other files.

In this case I am going to define the functions in such a way that it will take the control parameter `ControlFull`as one of the input, and then pass the corresponding values in the functions.

I thing eventually I will implement everything in one big function, but for now I will proceed step by step, this will be then used as a blueprint for future work.
=#
#=
#### 4.1. $ G_n $

To calculate the $ G_n $, I only need to calculate the integral I defined earlier and do that for different values of $n$.

I would like to define a general function that takes many arguments as parameter.
I know it is not the best way to do that, but it will help me reducing the allocation in the long term.
=#

function G_factor(n::Int64, z::Float64, η::ComplexF64, k::Float64, h::Float64)
    return time_dependent(n, η, k) *             # time dependent part
           spatial_fourier(n, η, z) *           # spatial part, of energy level n
           (bh_integrand(z, h, η) - sd_groundstate(z, η))
end

#=
#### 4.2. $ K_n $

The code to calculate the $ K_n $ is very similar to the one used to calculate the $ G_n $, the only big difference is that I have to repeat the calculation for the different terms of the gradient.

The nice part is that I do not have to use the function depending on the second derivative.

Similarly to what I have done with the $ G_n $ term, I will define a function that takes a lot of arguments.
=#

function K_factor(n::Int64, z::Float64, η::ComplexF64, k::Float64, h::Float64)
    return time_dependent(n, η, k) *             # time dependent part
           spatial_fourier(n, η, z) *           # spatial part, of energy level n
           (bh_integrand(z, h, η))
end

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
function corrections(k::Vector{ComplexF64}, g::ComplexF64)
    H = k * k' # Hessian matrix 
    v = g' * k |> real # Fidelity vector
    num = norm(v)^3 * v
    den = v' * H * v
    return num / den |> real
end
corrections(g::ComplexF64, k::Vector{ComplexF64}) = corrections(k, g)
"""
    corrections(n::Int64, c::Control; λs::Int64=5)
Return a tuple of the form '(Kₙ, Gₙ)', given the energy level `n` and the control parameter `c`.
Everything is defined in place as it is easier to implement.
"""
function corrections(ns::Vector{Int64}, c::Control; λs::Int64=5)
    b(t) = auxiliary(t, c)
    db(t) = ForwardDiff.derivative(b, t)
    ddb(t) = ForwardDiff.derivative(db, t)
    J = control_function(b, ddb, c)
    ξ0, U, h = EBJJ.scaling_ξ0(c), c.U, 1.0 / c.N # constants
    k(t::Float64) = EBJJ.interpolation_integral(c)(t)
    η(t::Float64) = (ξ0^2 / b(t)^2 - 2im * db(t) / (U * b(t)))
    gradient_functions = gradient(collect(0.0:c.T/(λs+1):c.T))
    grad(t::Float64) = [g(t) for g in gradient_functions]
    corr = zeros(Float64, λs)
    for i in 1:length(ns)
        global num = ns[i]
        # Function to evaluate the Kns
        k(v) = K_factor(num, v[2], η(v[1]), k(v[1]), h) * grad(v[1])
        # Function to evaluate the Gns
        g(v) = G_factor(num, v[2], η(v[1]), k(v[1]), h) * J(v[1])
        Kn::Vector{ComplexF64} = hcubature(k, [0.0, -1.0e1], [c.T, 1.0e1], atol=1.0e-5, maxevals=100000)[1]
        Gn::ComplexF64 = hcubature(g, [0.0, -1.0e1], [c.T, 1.0e1], atol=1.0e-5, maxevals=100000)[1]
        corr += corrections(Kn, Gn)
    end
    return corr
end
