using Test
using Qecsim.Model
using Qecsim:QecsimError
using Qecsim.BasicModels:BasicCode
using Qecsim.PauliTools:to_bsf

@testset "StabilizerCode" begin
    # 5-qubit code
    code = BasicCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"], ["XXXXX"], ["ZZZZZ"])
    @test logicals(code) == to_bsf(["XXXXX", "ZZZZZ"])
    @test validate(code) === nothing
    # Non-commuting stabilizers
    code = BasicCode(["XXXXI", "IXZZX", "XIXZZ", "ZXIXZ"], ["XXXXX"], ["ZZZZZ"])
    @test_throws QecsimError validate(code)
    # Non-commuting stabilizers with logicals
    code = BasicCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"], ["XXXII"], ["IIZZZ"])
    @test_throws QecsimError validate(code)
    # Commuting logicals
    code = BasicCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"], ["IIIII"], ["IIIII"])
    @test_throws QecsimError validate(code)
end
