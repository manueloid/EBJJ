#=
# Single Parameter testing

In the notes I already set up the calculations where only one parameter is used.

Just for reference, the formula of the Fourier transform of the STA wave function where the argument of the Gaussian part is set as the parameter $ \alpha $ is given by
$$
    \chi_n(z,t) = \left [ \frac{R_\alpha^2}{π}]^{1/4} \frac{1}{\sqrt{2^nn!} \exp\left\{- i \left(n + \frac{1}{2} \right) \int_0^t d\tau ~ U R_\alpha^2(\tau) \right\}
    \frac{i^n}{\sqrt{\hbar}}\frac{1}{\alpha} \left[\frac{(\alpha^2)^*}{\alpha^2}\right]^{n/2}
    \exp{-\frac{z^2}{2 \hbar^2 \alpha^2}}
    \mathcal{H}_n\left(\frac{R_\alpha}{|\alpha^2|} \frac{z}{\hbar}\right)
    \tag{1}
$$

where - in terms of the external bosonic Josephson Junction parameters - the value $ \alpha ^2 $ is defined as
$$
  \alpha^2 = \frac{ \sqrt{8 J_0 N / U}}{b^2} - \frac{2 i \dot{b}}{U b} 
  \tag{2}
$$

and the term $ R $ is defined as the square root of the real part of $ \alpha ^2 $, i.e.
$$
R_\alpha = \sqrt{\Re(\alpha^2)} = \left(\frac{8J_0 N}{U}\right)^{1/4} \frac{1}{b}
  \tag{3}
$$

What I want to test here is if the simplifications I made when dealing with the complex numbers are correct.
In particular, I first want to check if the following equality holds
$$
2 \frac{R_\alpha^2}{\alpha^2} - 1 = /frac{(\alpha^2)^*}{\alpha^2} 
\tag{4}
$$

Then I want to numerially check if the following relation is true
$$
    \frac{R_\alpha}{\alpha^2 \sqrt{2 \frac{R_\alpha^2}{\alpha^2} - 1}} = \frac{R_\alpha}{|\alpha^2|}
\tag{5}
$$

As I said I am going to do that numerically, first I will pass some random complex number to see if the relations are correct and only then I will pass the actual values of $ \alpha^2 $ for different times $t$.
=#
#=
## 1 - Testing the relations with random complex numbers

Here I am going to define some functions to test the relations I discussed earlier, but in this case the argument will be only a general complex number.

The plan is to test the relation for a certain number of general complex number and make sure that the relation holds.
=#

using Test
"""
    test_1(param::ComplexF64)
Given a complex number `param`, this function tests if the relation (4) holds.
The relation is 2 ℜ(α) / α - 1 = (α*) / α
"""
function test_1(param::ComplexF64)
    r = sqrt(real(param))
    result = 2 * r^2 / param - 1 - (conj(param) / param)
    return result
end
"""
    test_2(param::ComplexF64)
Given a complex number `param`, this function tests if the relation (5) holds.
The relation is α^2 * sqrt(2 * R_α^2 / α^2 - 1) = |α^2|
"""
function test_2(param::ComplexF64)
    result = (param * sqrt(conj(param)/param)) - abs(param)
    return result
end
@testset "testing relation 2" begin
    for _ in 1:10000
        param = rand(ComplexF64)
        @test isapprox(test_1(param), 0.0, atol=1e-10)
        @test isapprox(test_2(param), 0.0, atol=1e-10)
    end
end

#=
## 2 - Testing the relations with the actual values of α^2

Now it is time to pass the actual values of $ \alpha^2 $ and see if the relations hold.
There is no need to define any new functions, but only to get the value of $ \alpha^2 $ at a certain time $ t $ and then pass it to the functions I defined earlier.

The function α2 is the one that returns the value of $ \alpha^2 $ at a certain time $ t $.
This function will also accept some control parameter variable that are the ones that give the parameters of the system.

I would like to point out that in this case I did not make the substitution $ \hbar \rigtarrow 1 / N$.
=#

function α2(t::Float64, c::Control)
    J0, N, U = c.J0, c.N, c.U
    h = 1 / N
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    α(t::Float64) = ( sqrt(8J0 * N / U) / b(t)^2 - 2im / U * db(t) / b(t) ) 
    return α(t)
end
c = ControlFull()
t = range(0.0, c.T, length=10000)
@testset "testing relation 2" begin
    for i in t
        @test isapprox(test_1(α2(i, c)), 0.0, atol=1e-10)
        @test isapprox(test_2(α2(i, c)), 0.0, atol=1e-10)
    end
end

#=
## 3 - Takeaway

In the previous tests I have shown that the relations I used in the notes are correct.
Now I can write the main parts of the Fourier transform of the STA wave function using only one complex parameter.

In particular, I can write the whole wave function as the one in equation (1) at the top of this file.
=#
