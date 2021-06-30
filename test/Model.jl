using Test
using Qecsim.Model
using Qecsim.Models.Basic
using Qecsim:PauliTools as PT

@testset "StabilizerCode" begin
    # 5-qubit code
    code = BasicCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"], ["XXXXX"], ["ZZZZZ"])
    @test logicals(code) == PT.pauli_to_bsf(["XXXXX", "ZZZZZ"])
    @test validate(code) === nothing
    # Non-commuting stabilizers
    code = BasicCode(["XXXXI", "IXZZX", "XIXZZ", "ZXIXZ"], ["XXXXX"], ["ZZZZZ"])
    @test_throws AssertionError validate(code)
    # Non-commuting stabilizers with logicals
    code = BasicCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"], ["XXXII"], ["IIZZZ"])
    @test_throws AssertionError validate(code)
    # Commuting logicals
    code = BasicCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"], ["IIIII"], ["IIIII"])
    @test_throws AssertionError validate(code)
end
