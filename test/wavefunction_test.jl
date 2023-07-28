#=
# Wave function testing

In order to try to minimize the time it takes to calculate the corrections, I am going to interpolate the integral function and then define a new one.
I will define an array of points where I am going to sample the integral, and then I will interpolate the function at those points.
Then I will test if the values of the functions are similar up to a certain tolerance.
=#

using QuadGK
using Interpolations
using BenchmarkTools
using ForwardDiff
using SpecialPolynomials
using EBJJ

"""
    interpolation_check(c::Control; npoints=1000, checkpoints=10000, atol=1e-3 )
Return two lists of values

1. The first array contains the value of the interpolating function at each point in `range(0.0, c.T, length=checkpoints)`, given the interpolating function is defined on the grid `range(0.0, c.T, length=npoints)`
2. The second array contains the value of the integral function at each point in `range(0.0, c.T, length=checkpoints)`.
"""
function interpolation_check(c::Control; npoints=1000, checkpoints=10000)
    trange = range(0.0, c.T, length=checkpoints)
    itp = EBJJ.interpolation_integral(c, npoints=npoints)
    phase_integrand(t, c::Control) = 1 / auxiliary(t, c)^2
    int(t) = quadgk(t -> phase_integrand(t, c), 0, t)[1]
    interpolated_values = [itp(t) for t in trange]
    integral_values = [int(t) for t in trange]
    return interpolated_values, integral_values
end

c = ControlFull(10, 1.0)
interpolated_values, integral_values = interpolation_check(c, npoints=1000, checkpoints=10000)
@testset "Approximation of the integral function" begin
    @test isapprox(interpolated_values, integral_values, atol=1e-3)
end

#=
## 2. Thourough test of the wave functions.

The idea is to have the numerically implement the Fourier transform of the wave function and compare the result with the analytical solution.

### 2.1 Wave function with no simplification

The plan here is to focus on the single wave function without any simplification.
I will write the general STA wave function as the product of 5 terms:
1. the normalisation factor, which is the same as the one for the Harmonic Oscillator 
$$ \frac{\sqrt{\xi_0}}{\sqrt[4]{\pi }$$

2. the imaginary phase, given by 
$$ \exp \left(-i \left(n+\frac{1}{2}\right) \int_0^t \frac{\xi_0^2 U}{2 b(\tau )^2} \, d\tau \right) $$

3. the Gaussian term, given by
$$ \exp \left(\frac{1}{2} p^2 \left(-\frac{\xi_0^2}{b(t)^2}+\frac{2 i b'(t)}{U b(t)}\right)\right) $$

4. the Hermite polynomial, given by
$$ \mathcal{H}_n\left(\frac{p \xi_0}{b(t)}\right). $$ 

5. the factor depending on the excitiation, given by
$$\sqrt{2^n n! b(t)}}$$

The functions I can reuse from the main code are the following:
- `scaling_ξ0`
- `interpolation_integral`, which returns the interpolating function for the integral of the phase
- `imaginary_phase`, which returns the whole function of the imaginary phase given the index `n` and the time `t`

The functions I need to implement are the following:
- `gaussian`, which gives the value of the Gaussian term at a given time `t` and momentum `p`
- `hermite`, which gives the value of the Hermite polynomial at a given momentum `p`, index `n` and time `t`.
- `excitation`, which gives the value of the excitation term at a given index `n` and time `t`.

Then I am going to define the wave function as the product of all these terms.
At the same time I am going to define a big wave function function where all the terms are defined in the body, instead of being a product of different terms, to see if there is any difference.
Again, I will this is the momentum representation of the wave function, and what I want to check is that the normalisation is conserved.

=#
"""
    `normalisation(c::Control)`
Returns the normalisation factor of the wave function.
"""
normalisation(c::Control) = sqrt(EBJJ.scaling_ξ0(c)) / (π)^0.25
"""
    `gaussian(t, p, c::Control)`
Returns the value of the Gaussian term at a given time `t` and momentum `p`.
"""
function gaussian(t, p, c::Control)
    ξ0, U = EBJJ.scaling_ξ0(c), c.U # constants
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    return exp(0.5 * p^2 * (-ξ0^2 / b(t)^2 + 2im * db(t) / (U * b(t))))
end
"""
    `hermite(n, t, p c::Control)`
Returns the value of the Hermite polynomial at a given momentum `p`, index `n` and time `t`.
"""
function hermite(n, t, p, c::Control)
    ξ0 = EBJJ.scaling_ξ0(c) # constants
    b(t) = auxiliary(t, c) # Auxiliary function
    return SpecialPolynomials.basis(Hermite, n)(p * ξ0 / b(t)) # Hermite polynomial
end
"""
    `excitation(n, t, c::Control)`
Returns the value of the excitation term at a given index `n` and time `t`.
"""
function excitation(n, t, c::Control)
    b(t) = auxiliary(t, c) # Auxiliary function
    return sqrt(2^n * factorial(n) * b(t))
end
"""
    `wave_function(n, t, p, c::Control)`
This is the wave function in the momentum representation, with no simplification where all its parts are defined in its body
"""
function wave_function(n, t, p, c::Control)
    ξ0, U = EBJJ.scaling_ξ0(c), c.U # constants
    b(t) = auxiliary(t, c) # Auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # Derivative of the auxiliary function
    normalisation = sqrt(EBJJ.scaling_ξ0(c)) / (π)^0.25 # Normalisation factor
    itp = EBJJ.interpolation_integral(c) # Interpolating function for the integral of the phase
    imaginary_phase(n, t) = exp(-im * (n + 0.5) * ξ0^2 * U / 2 * itp(t)) # Imaginary phase
    gaussian(t, p) = exp(0.5 * p^2 * (-ξ0^2 / b(t)^2 + 2im * db(t) / (U * b(t)))) # Gaussian term
    hermite(n, t, p) = SpecialPolynomials.basis(Hermite, n)(p * ξ0 / b(t)) # Hermite polynomial
    excitation(n, t) = sqrt(2^n * factorial(n) * b(t)) # Excitation term
    return normalisation * imaginary_phase(n, t) * gaussian(t, p) * hermite(n, t, p) / excitation(n, t)
end
"""
    `wave_function_split(n, t, p, c::Control)`
This function is the product of the smaller functions I defined earlier.
"""
function wave_function_split(n, t, p, c::Control)
    ξ0, U = EBJJ.scaling_ξ0(c), c.U # constants
    itp = EBJJ.interpolation_integral(c) # Interpolating function for the integral of the phase
    imaginary_phase(n, t) = exp(-im * (n + 0.5) * ξ0^2 * U / 2 * itp(t)) # Imaginary phase
    return normalisation(c) * imaginary_phase(n, t) * gaussian(t, p, c) * hermite(n, t, p, c) / excitation(n, t, c)
end

@testset "Wave function test" begin
    c = ControlFull(10, 0.1)
    t = [0.0, c.T / 2, c.T]
    zrange = -10.0:0.1:10.0
    for n in 0:2
        for t in t
            for p in zrange
                @test wave_function(n, t, p, c) ≈ wave_function_split(n, t, p, c)
            end
        end
    end
end

#=
### 2.2 Benchmarking

In this part I am going to benchmark the two functions that I defined earlier.
I am going to use the `BenchmarkTools` package to do this.
I will integrate the two functions for a given time and then I will compare the time it takes to do so.
This will be done with the `QuadGK` package.
=#

c = ControlFull(10, 0.1)
t = c.T
# @btime quadgk(x -> wave_function(0, t, x, c), -10.0, 10.0)
# @btime quadgk(x -> wave_function_split(0, t, x, c), -10.0, 10.0)
# @btime quadgk(x -> wave_function(2, t, x, c), -10.0, 10.0)
# @btime quadgk(x -> wave_function_split(2, t, x, c), -10.0, 10.0)

#=
There is no clear difference in the time it takes to integrate the two functions.
But I think I will with the first one so that everything is self contained.

### 2.3 Normalisation
Now I need to check that the wave functions are normalised.
As I said, I will stick to the `wave_function` function, because it is more compact.
The normalisation will be done by integrating the square of the wave function over all space.

After that I will make a simplified version of the wave function I need to integrate, and check that it is normalised.
The simplified version of the wave function is the function $ |\chi_n(p, t)|^2 $, where we can see that the imaginary part cancels out and the Gaussian term is $e^{p^2}$, instead of $e^{p^2/2}$.

Finally I will benchmark the normalisation of the two functions, in theory the simplified version should be (a lot) faster.
=#
"""
    `normalisation(n, t, c::Control)`
Take the nth wave function, calculate the absolute value squared and integrate it over all space.
The result is given at a time `t`.
"""
function normalisation(n, t, c::Control)
    return quadgk(x -> conj(wave_function(n, t, x, c)) * wave_function(n, t, x, c), -10.0, 10.0)[1]
end
"""
    `wave_function_simplified(n, t, p, c::Control)`
Gives the value of the absolute value squared of the wave function, at a given time and momentum, with all the simplifications carried out.
"""
function wave_function_simplified(n, t, p, c::Control)
    ξ0 = EBJJ.scaling_ξ0(c) # constants (no U is needed in this case)
    b(t) = auxiliary(t, c) # Auxiliary function (no derivative is needed)
    normalisation = ξ0 / sqrt(π) # Normalisation factor
    gaussian(t, p) = exp(p^2 * (-ξ0^2 / b(t)^2)) # Gaussian term
    hermite(n, t, p) = SpecialPolynomials.basis(Hermite, n)(p * ξ0 / b(t)) # Hermite polynomial
    excitation(n, t) = 2^n * factorial(n) * b(t) # Excitation term
    return normalisation * gaussian(t, p) * hermite(n, t, p)^2 / excitation(n, t)
end
"""
    `normalisation_simplified(n, t, c::Control)`
Similar to `normalisation`, but using the simplified version of the wave function.
"""
function normalisation_simplified(n, t, c::Control)
    return quadgk(x -> wave_function_simplified(n, t, x, c), -1.0, 1.0)[1]
end
@testset "Normalisation momentum space" begin
    c = ControlFull(10, 0.1)
    t = [0.0, c.T / 2, c.T]
    n = 2
    for t in t
        for n in n
            @test isapprox(normalisation(n, t, c) |> real, 1.0)
            @test isapprox(normalisation_simplified(n, t, c), 1.0)
        end
    end
end
#=
### 2.3 Benchmarking
Now I can benchmark the two normalisation functions.
In theory the simplified version should be faster, because it doesn't have to calculate the imaginary part of the wave function.
=#
c = ControlFull(10, 0.1)
t = 0.0
n = 2
# @btime normalisation(n, t, c)
# @btime normalisation_simplified(2, t, c)
#=
As expected, the simplified version is faster.
To sum up, we had 
- `682.422 ms (214208 allocations: 53.29 MiB)` for the normalisation of the full wave function
- `61.522 μs (698 allocations: 55.77 KiB)` for the normalisation of the simplified wave function

### 2.4 Checking orthogonality
Just for a sanity check, I will check that the wave functions are orthogonal.
I will try only the full wave function, as I do not want ot implement the simplified version of the $\langle \chi_n | \chi_m \rangle$ integrand.
=#

"""
    `orthogonality_check(n,m, t, c::Control)`
Checks the orthogonality of the nth state and the mth state, by evaluating the integral from -1 to 1 by performing the integral over z 
"""
function orthogonality_check(n, m, t, c::Control)
    integrand(z) = conj(wave_function(n, t, z, c)) * wave_function(m, t, z, c)
    return quadgk(x -> integrand(x), -10.0, 10.0)[1] |> real
end

@testset "Orthogonality" begin
    c = ControlFull(10, 0.1)
    t = [0.0, c.T / 2, c.T]
    n = [0, 1, 2]
    m = [3, 4, 5]
    for t in t
        for (n, m) in zip(n, m)
            @test isapprox(orthogonality_check(n, m, t, c), 0.0)
        end
    end
end
