using Test
using Qecsim.BasicModels
using Qecsim.Model
using Qecsim.PauliTools:to_bsf


@testset "BasicCode" begin
    pauli_stabilizers = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]
    pauli_logical_xs = ["XXXXX"]
    pauli_logical_zs = ["ZZZZZ"]
    my_nkd = (5, 1, 3)
    my_label = "5-qubit"
    code = BasicCode(pauli_stabilizers, pauli_logical_xs, pauli_logical_zs, my_nkd,
                     my_label)
    @test nkd(code) == my_nkd
    @test label(code) == my_label
    @test stabilizers(code) == to_bsf(pauli_stabilizers)
    @test logical_xs(code) == to_bsf(pauli_logical_xs)
    @test logical_zs(code) == to_bsf(pauli_logical_zs)
    # defaults
    code = BasicCode(pauli_stabilizers, pauli_logical_xs, pauli_logical_zs)
    @test isequal(nkd(code), (5, 1, missing))  # isequal equates missing
    @test isa(label(code), String)
    # 5-qubit
    code = FiveQubitCode()
    @test nkd(code) == (5, 1, 3)
    @test label(code) == "5-qubit"
    @test validate(code) === nothing  # no error
    # Steane
    code = SteaneCode()
    @test nkd(code) == (7, 1, 3)
    @test label(code) == "Steane"
    @test validate(code) === nothing  # no error
end
