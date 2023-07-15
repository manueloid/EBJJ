using EBJJ
using Test

@testset "EBJJ.jl" begin
    @test ControlFull().N == 30
end
