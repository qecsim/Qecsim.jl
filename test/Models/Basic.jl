using Qecsim:Model
using Qecsim.Models:Basic

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
end
