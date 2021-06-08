using Qecsim: PauliTools as PT
using Test

@testset "PauliTools.jl" begin
    @test PT.hello("Boo") == "Hello Boo"

    @testset "bsp" begin
        # III bsp III commute
        @test PT.bsp(BitVector([0, 0, 0, 0, 0, 0])', BitVector([0, 0, 0, 0, 0, 0])) == 0
        # XII bsp XII commute
        @test PT.bsp(BitVector([1, 0, 0, 0, 0, 0])', BitVector([1, 0, 0, 0, 0, 0])) == 0
        # XII bsp IZZ commute
        @test PT.bsp(BitVector([1, 0, 0, 0, 0, 0])', BitVector([0, 0, 0, 0, 1, 1])) == 0
        # XII bsp ZII do not commute
        @test PT.bsp(BitVector([1, 0, 0, 0, 0, 0])', BitVector([0, 0, 0, 1, 0, 0])) == 1
        # XXI bsp ZZI commute
        @test PT.bsp(BitVector([1, 1, 0, 0, 0, 0])', BitVector([1, 1, 0, 0, 0, 0])) == 0
        # XXX bsp ZZZ do not commute
        @test PT.bsp(BitVector([1, 1, 1, 0, 0, 0])', BitVector([0, 0, 0, 1, 1, 1])) == 1

        # XII bsp IYY commute
        @test PT.bsp(BitVector([1, 0, 0, 0, 0, 0])', BitVector([0, 1, 1, 0, 1, 1])) == 0
        # XII bsp YII do not commute
        @test PT.bsp(BitVector([1, 0, 0, 0, 0, 0])', BitVector([1, 0, 0, 1, 0, 0])) == 1
        # XXI bsp YYI commute
        @test PT.bsp(BitVector([1, 1, 0, 0, 0, 0])', BitVector([1, 1, 0, 1, 1, 0])) == 0
        # XXX bsp YYY do not commute
        @test PT.bsp(BitVector([1, 1, 1, 0, 0, 0])', BitVector([1, 1, 1, 1, 1, 1])) == 1

        stabilizers = BitMatrix(  # 5-qubit stabilizers
            [1 0 0 1 0 0 1 1 0 0    # XZZXI
             0 1 0 0 1 0 0 1 1 0    # IXZZX
             1 0 1 0 0 0 0 0 1 1    # XIXZZ
             0 1 0 1 0 1 0 0 0 1])  # ZXIXZ
        error = BitVector([0, 0, 1, 1, 0, 0, 1, 0, 1, 0])  # IZXYI
        commutation = [0, 1, 1, 0]
        @test PT.bsp(stabilizers, error) == commutation
        @test PT.bsp(error', stabilizers') == commutation'
        errors = hcat(
            BitVector([1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),  # XIIII
            BitVector([0, 0, 0, 0, 0, 0, 0, 1, 0, 0]),  # IIZII
            BitVector([0, 0, 0, 0, 1, 0, 0, 0, 0, 1]),  # IIIIY
            BitVector([0, 0, 1, 1, 0, 0, 1, 0, 1, 0]),  # IZXYI
        )
        commutations = [[0, 0, 0, 1] [0, 0, 1, 0] [0, 1, 1, 1] [0, 1, 1, 0]]
        @test PT.bsp(stabilizers, errors) == commutations
        @test PT.bsp(errors', stabilizers') == commutations'
    end
end
