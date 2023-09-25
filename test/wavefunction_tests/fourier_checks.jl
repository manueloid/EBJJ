#=
# Fourier checks

In this file I am going to numerically check if the analytic form of the Fourier transform of the STA wave function is correct.
The plan is to numerically implement the Fourier transform and then pass the same value to the analytic form of the Fourier transform and check if the two are equal.

I will first do that for a general complex number and only later I will pass the actual value of the parameter $ \alpha^2$.

I already tested that the analytic solution can be expressed in terms of only on parameter $ \alpha^2$ and its real part $ R_\alpha$.
I will use this fact to simplify the code.
=#
#=
### 1 Fourier transform general complex number

The function I need to find the Fourier transform of is the following:
$$
    \exp\left\{	- \frac{p^2}{2} \alpha^2\right\}
    \mathcal{H}_n\left( p R_/alpha \right)
    \tag{1}
$$

where in this case the parameter $ \alpha $ is not dependent on time $ t $, this will be discussed later.
In the context of quantum mechanics, the Fourier transform is quite different as we need an extra term $ \hbar$ in the exponent and in the normalisation term.
In particular, the integral is
$$ 
	\frac{1}{\sqrt{2\pi\hbar}}
	\int_{-\infty}^{+\infty}
    \exp\left\{	- \frac{p^2}{2} \alpha^2\right\}
    \mathcal{H}_n\left( p R_/alpha \right)
	e^{\frac{i}{\hbar}px} dp
    \tag{2}
$$

And we can simplify this expression by defining a new variable  $ \tilde{z} = z / \hbar$.

The analytic solution of this integral is given by
$$
    \frac{i^n}{\sqrt{\hbar}} 
    \frac{1}{\alpha}
    \left[\frac{\alpha^{2*}}{\alpha^2}\right]^{n/2}
    \exp\left\{- \frac{\tilde{z}^2}{2\alpha^2}\right\}
    \mathcal{H}_n\left( \tilde{z} \frac{R_\alpha}{|\alpha^2} \right)
    \tag{3}
$$

I am going to test this relation for the new defined variable $ \tilde{z} $, I can always pass the actual value of $ z $ later.

#### 1.1 Code implementation

In terms of the implementation of the code, I will define two versions of each term of the function.
One will be the numerical solution to the integral and I will add the suffix `_num`, while the other one will be the analytical solution of the integral and I will use the suffix `_an`.

I will only care about the Gaussian term and the Hermite polynomial, as the other terms do not depend on the variable $ z $.
Moreover, I will only define the full function, I will not split the whole function in smaller chunks.
This is because this file is only a test and I do not need any optimisation.
=#

he(n::Int64, z::Float64) = SpecialPolynomials.basis(Hermite, n)(z)
"""
    toint(n::Int64, p::Float64, α2::ComplexF64)
This function is the integrand of the Fourier transform of the STA wave function.
It is basically only the product of the Gaussian term and the Hermite polynomial.
The function return the value of the integrand for the momentum `p` and for the parameter `α2`, given the excitation level `n`.
Is the one I am going to pass to the `quadgk` function.
"""
function toint(n::Int64, p::Float64, α2::ComplexF64)
    r = sqrt(real(α2))
    return exp(-p^2 * α2 / 2) * he(n, p * r)
end
"""
    fourier(n::Int64, z::Float64, α2::ComplexF64)
Numerically return the value of the Fourier transform of `toint` for the excitation level `n`, the position `z` and the parameter `α2`.
It also takes the value of the parameter `h` which is the Planck constant.
"""
function fourier_num(n::Int64, z::Float64, α2::ComplexF64, h::Float64)
    return quadgk(p -> toint(n, p, α2) * exp(im * p * z / h), -Inf, Inf, atol=1e-7)[1] / sqrt(2 * π * h)
end
"""
    fourier_an(n::Int64, z::Float64, α2::ComplexF64, h::Float64)
Return the value of the analytic solution of the Fourier transform of the STA wave function for the excitation level `n`, the position `z` and the parameter `α2`.
"""
function fourier_an(n::Int64, z::Float64, α2::ComplexF64, h::Float64)
    r = sqrt(real(α2))
    return im^n / sqrt(h) * 1 / sqrt(α2) * (conj(α2) / α2)^(n / 2) * exp(-z^2 / (2 * α2)) * he(n, z * r / abs(α2))
end
@testset "Fourier transform general complex" begin
    for _ in 1:10000
        n = rand(0:4)
        z = rand()
        α2 = rand(ComplexF64)
        h = 0.01
        @test isapprox(fourier_num(n, z, α2, h), fourier_an(n, z / h, α2, h), atol=1e-4)
    end
end

#=
### 2 Fourier transform actual values 

In the next section I am going to pass the actual value of the parameter $ \alpha^2 $, as I have already done in other files.
=#

function α2(t::Float64, c::Control)
    J0, N, U = c.J0, c.N, c.U
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    α(t::Float64) = (sqrt(8J0 * N / U) / b(t)^2 - 2im / U * db(t) / b(t))
    return α(t)
end
@testset "Fourier transform actual values" begin
    for _ in 1:1000
        n = rand(0:4)
        N = rand(10:10:50)
        c = ControlFull(N, 0.03)
        h = 1 / N
        t = rand(0.0:c.T)
        @test isapprox(fourier_num(n, z, α2(t, c), h), fourier_an(n, z / h, α2(t, c), h), atol=1e-4)
    end
end
