using Qecsim
using Test

@testset "Qecsim.jl" begin
    @test Qecsim.doubler(3) == 6
end
