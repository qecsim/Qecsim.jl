using Qecsim:Model
using Qecsim.Models:Basic
using Qecsim:PauliTools as PT


@testset "BasicCode" begin
    pauli_stabilizers = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]
    pauli_logical_xs = ["XXXXX"]
    pauli_logical_zs = ["ZZZZZ"]
    nkd = (5, 1, 3)
    label = "5-qubit"
    code = Basic.BasicCode(pauli_stabilizers, pauli_logical_xs, pauli_logical_zs;
        nkd=nkd, label=label)
    @test Model.nkd(code) == nkd
    @test Model.label(code) == label
    @test Model.stabilizers(code) == PT.pauli_to_bsf(pauli_stabilizers)
    @test Model.logical_xs(code) == PT.pauli_to_bsf(pauli_logical_xs)
    @test Model.logical_zs(code) == PT.pauli_to_bsf(pauli_logical_zs)
    # defaults
    code = Basic.BasicCode(pauli_stabilizers, pauli_logical_xs, pauli_logical_zs)
    @test Model.nkd(code) == (5, 1, nothing)
    @test isa(Model.label(code), String)
    # 5-qubit
    code = Basic.FiveQubitCode()
    @test Model.nkd(code) == (5, 1, 3)
    @test Model.label(code) == "5-qubit"
    @test Model.validate(code) === nothing  # no error
    # Steane
    code = Basic.SteaneCode()
    @test Model.nkd(code) == (7, 1, 3)
    @test Model.label(code) == "Steane"
    @test Model.validate(code) === nothing  # no error
end
