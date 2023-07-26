#=
# Test suite for checking the control functions in this package.

=#
@testset "Polynomial test" begin
    c = ControlFull()
    @test control_function(-0.25, c) ≈ c.J0 # check if p(t) for t < 0 is J0
    @test control_function(0.0, c) ≈ c.J0 # check if p(0) is J0
    @test control_function(c.T, c) ≈ c.Jf # check if p(tf) is Jf
    @test control_function(1.25 * c.T, c) ≈ c.Jf # check if p(t) for t > tf is Jf
end
