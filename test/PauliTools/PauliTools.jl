import Qecsim.PauliTools as PT
using Test

@testset "PauliTools.jl" begin
    @test PT.hello("Boo") == "Hello Boo"
end
