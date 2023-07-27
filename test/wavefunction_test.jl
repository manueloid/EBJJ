#=
# Wave function testing

In order to try to minimize the time it takes to calculate the corrections, I am going to interpolate the integral function and then define a new one.
I will define an array of points where I am going to sample the integral, and then I will interpolate the function at those points.
Then I will test if the values of the functions are similar up to a certain tolerance.
=#

using QuadGK
using Interpolations
using BenchmarkTools
using EBJJ

"""
    interpolation_check(c::Control; npoints=1000, checkpoints=10000, atol=1e-3 )
Return two lists of values

1. The first array contains the value of the interpolating function at each point in `range(0.0, c.T, length=checkpoints)`, given the interpolating function is defined on the grid `range(0.0, c.T, length=npoints)`
2. The second array contains the value of the integral function at each point in `range(0.0, c.T, length=checkpoints)`.
"""
function interpolation_check(c::Control; npoints=1000, checkpoints=10000)
    trange = range(0.0, c.T, length=10000)
    itp = EBJJ.interpolation_integral(c, npoints=1000)
    phase_integrand(t, c::Control) = 1 / auxiliary(t, c)^2
    int(t) = quadgk(t -> phase_integrand(t, c), 0, t)[1]
    interpolated_values = [itp(t) for t in trange]
    integral_values = [int(t) for t in trange]
    return interpolated_values, integral_values
end

c = ControlFull(10, 1.0)
interpolated_values, integral_values = interpolation_check(c, npoints=1000, checkpoints=10000)
@test isapprox(interpolated_values, integral_values, atol=1e-3)

#=
## 2. Checking orthogonality

The functions that I implemented, make it extremely easy to check the orthogonality of the wave functions.
In particular, I can just fix a time and then take the integral of the two spatial parts of the wave functions.
This should be 0, I shouldn't even need to use the purely time-dependent terms.
=#

"""
    `orthogonality_check(n, t, c::Control)`
Checks the orthogonality of the nth state and the ground state, by evaluating the integral from -1 to 1.
"""
function orthogonality_check(n, t, c::Control)
    integrand(x) = conj(ground_state(t, x, c)) * ground_state(t, x, c)
    return quadgk(x -> integrand(x), -10.0, 10.0)[1]
end

@testset "Orthogonality" begin
    c = ControlFull(10, 0.1)
    @test isapprox(orthogonality_check(1, 0.0, c), 0.0, atol=1e-3)
    @test isapprox(orthogonality_check(2, 0.0, c), 0.0, atol=1e-3)
    @test isapprox(orthogonality_check(3, 0.0, c), 0.0, atol=1e-3)
    @test isapprox(orthogonality_check(4, 0.0, c), 0.0, atol=1e-3)
end

test
