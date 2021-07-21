using Test

using Qecsim.GenericModels
using Qecsim.Model
using Qecsim.BasicModels:FiveQubitCode, SteaneCode
using Qecsim.PauliTools:bsp, to_bsf

@testset "NaiveDecoder" begin
    decoder = NaiveDecoder()
    # label
    @test label(decoder) == "Naive"
    # decode
    code = FiveQubitCode()
    for error_pauli in ("IIIII", "IXIII", "IIYII", "IIIIZ", "IZYXI", "YZYXX")
        error = to_bsf(error_pauli)
        syndrome = bsp(stabilizers(code), error)
        result = decode(decoder, code, syndrome)
        recovery = result.recovery
        @test bsp(stabilizers(code), recovery) == syndrome
        @test !any(bsp(stabilizers(code), xor.(error, recovery)))
    end
    # max_qubits
    code = SteaneCode()
    error = to_bsf("IIIIIII")
    syndrome = bsp(stabilizers(code), error)
    @test_throws ArgumentError decode(NaiveDecoder(5), code, syndrome)  # max_qubits=5
    @test decode(NaiveDecoder(7), code, syndrome).recovery == error  # max_qubits=7
    @test decode(NaiveDecoder(-1), code, syndrome).recovery == error  # max_qubits=unlimited
end
