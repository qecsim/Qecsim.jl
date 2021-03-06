using Test
using Qecsim.PauliTools

@testset "bsp" begin
    # III bsp III commute (i.e. return false)
    @test bsp(BitVector([0, 0, 0, 0, 0, 0])', BitVector([0, 0, 0, 0, 0, 0]))::Bool == 0
    # XII bsp XII commute
    @test bsp(BitVector([1, 0, 0, 0, 0, 0])', BitVector([1, 0, 0, 0, 0, 0]))::Bool == 0
    # XII bsp IZZ commute
    @test bsp(BitVector([1, 0, 0, 0, 0, 0])', BitVector([0, 0, 0, 0, 1, 1]))::Bool == 0
    # XII bsp ZII do not commute (i.e. return true)
    @test bsp(BitVector([1, 0, 0, 0, 0, 0])', BitVector([0, 0, 0, 1, 0, 0]))::Bool == 1
    # XXI bsp ZZI commute
    @test bsp(BitVector([1, 1, 0, 0, 0, 0])', BitVector([1, 1, 0, 0, 0, 0]))::Bool == 0
    # XXX bsp ZZZ do not commute
    @test bsp(BitVector([1, 1, 1, 0, 0, 0])', BitVector([0, 0, 0, 1, 1, 1]))::Bool == 1

    # XII bsp IYY commute
    @test bsp(BitVector([1, 0, 0, 0, 0, 0])', BitVector([0, 1, 1, 0, 1, 1]))::Bool == 0
    # XII bsp YII do not commute
    @test bsp(BitVector([1, 0, 0, 0, 0, 0])', BitVector([1, 0, 0, 1, 0, 0]))::Bool == 1
    # XXI bsp YYI commute
    @test bsp(BitVector([1, 1, 0, 0, 0, 0])', BitVector([1, 1, 0, 1, 1, 0]))::Bool == 0
    # XXX bsp YYY do not commute
    @test bsp(BitVector([1, 1, 1, 0, 0, 0])', BitVector([1, 1, 1, 1, 1, 1]))::Bool == 1

    stabilizers = BitMatrix(  # 5-qubit stabilizers
        [1 0 0 1 0 0 1 1 0 0   # XZZXI
         0 1 0 0 1 0 0 1 1 0   # IXZZX
         1 0 1 0 0 0 0 0 1 1   # XIXZZ
         0 1 0 1 0 1 0 0 0 1]  # ZXIXZ
    )
    error = BitVector([0, 0, 1, 1, 0, 0, 1, 0, 1, 0])  # IZXYI
    commutation = BitVector([0, 1, 1, 0])
    @test bsp(stabilizers, error)::BitVector == commutation
    @test bsp(error', stabilizers')::BitMatrix == commutation'
    errors = hcat(
        BitVector([1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),  # XIIII
        BitVector([0, 0, 0, 0, 0, 0, 0, 1, 0, 0]),  # IIZII
        BitVector([0, 0, 0, 0, 1, 0, 0, 0, 0, 1]),  # IIIIY
        BitVector([0, 0, 1, 1, 0, 0, 1, 0, 1, 0]),  # IZXYI
    )
    commutation = BitMatrix([[0, 0, 0, 1] [0, 0, 1, 0] [0, 1, 1, 1] [0, 1, 1, 0]])
    @test bsp(stabilizers, errors)::BitMatrix == commutation
    @test bsp(errors', stabilizers')::BitMatrix == commutation'
end

@testset "pack_unpack" begin
    # specific test
    bsf = to_bsf("XIZIY")  # 1000100101
    packed_bsf = pack(bsf)
    @test packed_bsf == ("225", 10)
    @test unpack(packed_bsf) == bsf
    # long (>128) bit-vector (to test BigInt handles long binary vectors)
    bsf = BitVector(rand(0:1, 150))
    @test unpack(pack(bsf)) == bsf
end

@testset "to_bsf" begin
    # Single Paulis
    @test to_bsf("XIZIY") == BitVector([1, 0, 0, 0, 1, 0, 0, 1, 0, 1])
    @test to_bsf("IIIII") == BitVector([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    @test to_bsf("XXXXX") == BitVector([1, 1, 1, 1, 1, 0, 0, 0, 0, 0])
    @test to_bsf("ZZZZZ") == BitVector([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
    @test to_bsf("YYYYY") == BitVector([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
    # Multiple Paulis
    @test to_bsf(["XIZIY", "IXZYI"]) == BitMatrix(
        [1 0 0 0 1 0 0 1 0 1
         0 1 0 1 0 0 0 1 1 0]
    )
    @test to_bsf(["XXXXX", "ZZZZZ", "YYYYY"]) == BitMatrix(
        [1 1 1 1 1 0 0 0 0 0
         0 0 0 0 0 1 1 1 1 1
         1 1 1 1 1 1 1 1 1 1]
    )
end

@testset "to_pauli" begin
    # Single bsf
    @test to_pauli(BitVector([1, 0, 0, 0, 1, 0, 0, 1, 0, 1])) == "XIZIY"
    @test to_pauli(BitVector([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])) == "IIIII"
    @test to_pauli(BitVector([1, 1, 1, 1, 1, 0, 0, 0, 0, 0])) == "XXXXX"
    @test to_pauli(BitVector([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])) == "ZZZZZ"
    @test to_pauli(BitVector([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])) == "YYYYY"
    # Multiple bsfs
    @test to_pauli(BitMatrix(
        [1 0 0 0 1 0 0 1 0 1
         0 1 0 1 0 0 0 1 1 0])
    ) == ["XIZIY", "IXZYI"]
    @test to_pauli(BitMatrix(
        [1 1 1 1 1 0 0 0 0 0
         0 0 0 0 0 1 1 1 1 1
         1 1 1 1 1 1 1 1 1 1])
    ) == ["XXXXX", "ZZZZZ", "YYYYY"]
    end

@testset "weight" begin
    for (pauli, wt) in (
        ("IIIII", 0),
        ("XIZIY", 3),
        ("XXXXX", 5),
        ("YYYYY", 5),
        ("ZZZZZ", 5),
        (["XIZIY", "XXXXX"], 8),
    )
        @test weight(pauli) == wt
        @test weight(to_bsf(pauli)) == wt
    end
end
