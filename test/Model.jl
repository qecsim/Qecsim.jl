using Qecsim:Model
using Qecsim.Models:Basic
using Qecsim:PauliTools as PT

@testset "StabilizerCode" begin
    # 5-qubit code
    code = Basic.BasicCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"], ["XXXXX"], ["ZZZZZ"])
    @test Model.logicals(code) == PT.pauli_to_bsf(["XXXXX", "ZZZZZ"])
    @test Model.validate(code) === nothing
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
