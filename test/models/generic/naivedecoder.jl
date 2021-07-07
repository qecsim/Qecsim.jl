using Test

using Qecsim.GenericModels
using Qecsim.Model
using Qecsim.BasicModels:FiveQubitCode
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
        @test all(bsp(stabilizers(code), xor.(error, recovery)) .== 0)
    end
end
