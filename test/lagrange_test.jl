#=
# Testing the Lagrange polynomials

I need to calculate the $ i $ -th term of the gradient, in order to evaluate the term $ K_n $.
This term is nothing more than the $ i $ -th Lagrange polynomial term over the interval $ [0, t_f] $, where the value at that point is $ 1 $ and all the other values are $ 0 $.

In general, a Lagrange polynomial is defined as the one that takes the values `[y₁, y₂, ..., yₙ]` at the points `[x₁, x₂, ..., xₙ]`.
In my case then, the interval `[x₁, x₂, ..., xₙ]` is the interval `[0, ..., tf]`, where the length of the second array is `n`.
The only boundary conditions I have here is that the value of the polynomial for the initial and final points are equal to 0.

I then have two options then:

1. Define a Lagrange polynomial that will assume the values $ 1 $ at each point of the array `[0, ..., tf]` and extract the $ i $ -th term of the polynomial.

2. Define a Lagrange polynomial that will assume the values $ 1 $ at the point `i` and $ 0 $ at all the other points and do that for each point of the array `[0, ..., tf]`.

I will try to do both and see which one is faster.

The first option is the hardest to implement as it looks like is not so easy to extract the $ i $ -th term of the polynomial.
The second one is quite simpler, the problem is that I have to define the polynomial for each point of the array, which is not so efficient.

Overall though, the second option looks like the best one, as it is extremely easy to implement.
=#
#=
### 1. Breakdown of the function

#### 1.1. Lagrange polynomial theory

In the eSTA scheme, the correction of the control function is given by a polynomial $ P(t)$ of $ \nu $ parameters $ \lambda_i $.
The boundary conditions for this correction polynomial is given by $ P(0) = P(t_f) = 0 $.
If we want to define the correction polynomial in the Lagrange basis, we can define it as 
$$
    P(t) = \sum_{i=1}^{\nu} \lambda_i \prod_{\substack{k=1 \\ k\neq i}}^{\nu} \frac{t - t_k}{t_i - t_k}
    \tag{1}
$$

where $ t_i $ are the points where the polynomial is evaluated.
In this case, $ t_i $ are equally spaced points in the interval $ [0, t_f] $.

The derivative of this polynomial with respect to the $ \lambda_i $ terms then is nothing more than the $i$ -th basis of the  Lagrange polynomial evaluated at the point $ t_i $.
This corresponds to have a Lagrange polynomial that assumes the value $ 1 $ at the point $ t_i $ and $ 0 $ at all the other points.
This process of course can be skipped for the initial and final points of the interval, as the polynomial will be equal to zero at those points but regardless, we need to include these points in the calculation and ensure the polynomial is 0.


#### 1.2. Implementation

The number of coefficients $ \lambda_i $ is not equal to the number of points in the interval $ [0, t_f] $, as the polynomial is zero at the initial and final points.
So, if we define `ncoeffs` as the number of coefficients we want our polynomial to have, the number of points in the interval will be `ncoeffs + 2`.

I think the best and most flexible way to implement this is to pass the array of points in the interval as an argument of the function, and retrieve the number of coefficients from it.

So, given an array of the form `[0, …, tf]`, the number of coefficients will be `length(array) - 2`.
The rest of the algorithm is shown in the definition of the function below.

Finally, I think the best output is an array of anonymous functions, where each function is the $ i $ -th term of the gradient.

=#

using Polynomials, SpecialPolynomials

#=
### 2. Testing the function

If everything goes according to plan, the functions evaluated at time `t` should return `1.0` for the `i` -th term and `0.0` for all the other terms.
Hence if I am going to loop over the array of functions and evaluate them for each time in the time array `tarray`, I should get the identity matrix as the output.

I am going to do that for an insane amount of points, just to be sure that the function is working properly.
I will also use the control parameter `ControlFull()` to mimic the calculations I will do later on.

I am not going to focus on the speed of the function, as it is only a test function and I do not need it to be efficient.
=#

using EBJJ, LinearAlgebra, Test

@testset "Testing the gradient function" begin
    c = ControlFull()
    xarr = range(0.0, stop=c.T, length=100) |> collect
    sol = EBJJ.gradient(xarr)
    id = Matrix(1.0I, length(xarr) - 2, length(xarr) - 2)
    @test [solution(x) for solution in sol, x in xarr[2:end-1]] ≈ id
end
