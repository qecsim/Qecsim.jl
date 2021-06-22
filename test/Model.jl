using Qecsim:Model
using Qecsim.Models:Basic

@testset "StabilizerCode" begin
    label = "my-basic-code"
    code = Basic.BasicCode(label)
    @test code.label == label
    @test Model.label(code) == label
    @test Model.label(Basic.BasicCode()) == "BasicCode"
end
