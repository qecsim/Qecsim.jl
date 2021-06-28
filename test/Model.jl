using Qecsim:Model
using Qecsim.Models:Basic

@testset "StabilizerCode" begin
    # 5-qubit code
    code = Basic.BasicCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"], ["XXXXX"], ["ZZZZZ"])
    Model.validate(code)
    # Non-commuting stabilizers
    code = Basic.BasicCode(["XXXXI", "IXZZX", "XIXZZ", "ZXIXZ"], ["XXXXX"], ["ZZZZZ"])
    @test_throws AssertionError Model.validate(code)
    # Non-commuting stabilizers with logicals
    code = Basic.BasicCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"], ["XXXII"], ["IIZZZ"])
    @test_throws AssertionError Model.validate(code)
    # Commuting logicals
    code = Basic.BasicCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"], ["IIIII"], ["IIIII"])
    @test_throws AssertionError Model.validate(code)
end
