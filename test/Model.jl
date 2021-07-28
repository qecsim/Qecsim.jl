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

@testset "DecodeResult" begin
    # typical valid kwargs/args constructors no exception
    # recovery only (typical)
    DecodeResult(recovery=BitVector([0, 0, 0, 0]))
    DecodeResult(nothing, BitVector([0, 0, 0, 0]), nothing, nothing)
    # success only (override)
    DecodeResult(success=true)
    DecodeResult(true, nothing, nothing, nothing)
    # success + logical_commutations (override)
    DecodeResult(success=true, logical_commutations=BitVector([0, 0]))
    DecodeResult(true, nothing, BitVector([0, 0]), nothing)
    # unusual but valid kwargs/args constructors no exception, type inferred
    # success + recovery (unusual)
    DecodeResult(success=true, recovery=[0, 0, 0, 0])
    DecodeResult(true, [0, 0, 0, 0], nothing, nothing)
    # recovery + logical_commutations (unusual)
    DecodeResult(recovery=[0, 0, 0, 0], logical_commutations=[0, 0])
    DecodeResult(nothing, [0, 0, 0, 0], [0, 0], nothing)
    # success + recovery + logical_commutations (unusual)
    DecodeResult(success=true, recovery=[0, 0, 0, 0],
        logical_commutations=[0, 0])
    DecodeResult(true, [0, 0, 0, 0], [0, 0], nothing)
    # custom_values valid kwargs/args constructors test type inferred
    # unspecified custom_values
    @test isa(DecodeResult(success=true), DecodeResult{Nothing})
    @test isa(DecodeResult(true, nothing, nothing, nothing), DecodeResult{Nothing})
    # custom_values Int
    @test isa(DecodeResult(success=true, custom_values=[1, 2]), DecodeResult{Vector{Int}})
    @test isa(DecodeResult(true, nothing, nothing, [1, 2]), DecodeResult{Vector{Int}})
    # custom_values Float64
    @test isa(DecodeResult(success=true, custom_values=[1., 2.]),
        DecodeResult{Vector{Float64}})
    @test isa(DecodeResult(true, nothing, nothing, [1., 2.]), DecodeResult{Vector{Float64}})
    # custom_values Int + Float64 promoted
    @test isa(DecodeResult(success=true, custom_values=[1, 2.]),
        DecodeResult{Vector{Float64}})
    @test isa(DecodeResult(true, nothing, nothing, [1, 2.]), DecodeResult{Vector{Float64}})
    # invalid kwargs/args constructors throw exception
    @test_throws QecsimError DecodeResult()
    @test_throws QecsimError DecodeResult(nothing, nothing, nothing, nothing)
    @test_throws QecsimError DecodeResult(logical_commutations=[1, 0])
    @test_throws QecsimError DecodeResult(nothing, nothing, [1, 0], nothing)
end
