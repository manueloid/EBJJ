#=
# Testing the corrections

I need to find out where the errors are when evaluating the corrections.

The problem I see is that the corrections do not scale well when the final time increases.
I already checked the auxiliary function and the problem is not there.

I will try to use two different integration methods, to see if the problem lies in there.

I think that generally I will implement functions that return the value of each term $ G_n $ and $ K_n$, so that I can use them in the future.

Moreover, since I do not need any optimisation, I will use the full definition of the STA wave function in position representation.
=#
#=
## 1 - Wave function 
=#

"""
    wave(z::Float64, t::Float64, c::Control)
Return the value of the STA wave function in position representation, for a given value of `z` and `t`.
It internally evaluates all the terms.
"""
function wave(z::Float64, t::Float64, c::Control)
    J0, Jf, U, tf, N = c.J0, c.Jf, c.U, c.tf, c.N
    b(t) = auxiliary(t, c)
    db(t) = auxiliary_1d(t, c)
    α2(t::Float64) = (1 / b(t)^2 * sqrt(8J0 * N / U) - 2im * db(t) / (U * b(t)))^(-1)  # Parameter of the Gaussian term
    φ(t::Float64) = √(2 * J0 * U * N)
    return
    (sqrt(8J0 * N / U) / π)^(1 / 4) *   #  normalisation term

end


#=
## 2 - Implementation using the `QuadGK` package

I will integrate over the two variables separately, first over `z` and then over `t`.
As `QuadGK` can integrate over infinite intervals, I will use that to my advantage.
I will check what is the numerical difference between integrating over `z ∈ [-Inf, Inf]` and `z ∈ [-1, 1]` and then over `z ∈ [-10, 10]`.
=#

#=
## 3 - Implementation using the `HCubature` package

Here I will not split the integration into two parts, but I will integrate over the variables `z` and `t` at the same time.
The main difference between this approach and the `QuadGK` one is that `HCubature` cannot integrate over infinite intervals, so I need to use an actual value.

I found that the best value for the moment is `z ∈ [-10, 10]` as it is the value that gives normalisation one for the `z` integral.
=#


