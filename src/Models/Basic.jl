"""
Basic stabilizer codes
"""
module Basic

import Qecsim.Model

using Qecsim.Model:StabilizerCode
using Qecsim:PauliTools as PT

export BasicCode, FiveQubitCode, SteaneCode

"""
    BasicCode <: StabilizerCode

Basic code defined by its stabilizers and logical operators.
"""
struct BasicCode <: StabilizerCode
    stablizers::BitMatrix
    logical_xs::BitMatrix
    logical_zs::BitMatrix
    nkd::Tuple{Int, Int, Union{Int,Nothing}}
    label::String
end
"""
    BasicCode(pauli_stabilizers, pauli_logical_xs, pauli_logical_zs; nkd=nothing,
        label=nothing)

Construct a basic code from string representations of stabilizers and logical operators.

Paulis are expressed as strings of capitalized I, X, Y, Z characters, with one character per
physical qubit. Logical X and Z operators are in matching order, with one of each for each
logical qubit. Optional `nkd` defaults to `n` and `k` evaluated and `d` nothing. Optional
`label` defaults to "Basic [n, k, d]".

# Examples
```jldoctest
julia> using Qecsim.Models.Basic

julia> code = BasicCode(["ZZI", "IZZ"], ["XXX"], ["IIZ"]);  # 3-qubit repetition

julia> validate(code)  # no error indicates operators satisfy commutation relations

julia> nkd(code)  # default nkd
(3, 1, nothing)

julia> label(code)  # default label
"Basic [3,1,nothing]"
```
"""
function BasicCode(pauli_stabilizers, pauli_logical_xs, pauli_logical_zs;
                   nkd=nothing, label=nothing)
    nkd = !isnothing(nkd) ? nkd : (
        length(pauli_stabilizers) > 0 ? length(pauli_stabilizers[1]) : 0,
        length(pauli_logical_xs),
        nothing)
    label = !isnothing(label) ? label : "Basic [$(nkd[1]),$(nkd[2]),$(nkd[3])]"
    BasicCode(PT.pauli_to_bsf(pauli_stabilizers), PT.pauli_to_bsf(pauli_logical_xs),
              PT.pauli_to_bsf(pauli_logical_zs), nkd, label)
end
Model.stabilizers(code::BasicCode) = code.stablizers
Model.logical_xs(code::BasicCode) = code.logical_xs
Model.logical_zs(code::BasicCode) = code.logical_zs
Model.nkd(code::BasicCode) = code.nkd
Model.label(code::BasicCode) = code.label

"""
    FiveQubitCode()

Construct 5-qubit [5,1,3] code as a [`BasicCode`](@ref).
"""
function FiveQubitCode()
    BasicCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"],
        ["XXXXX"], ["ZZZZZ"]; nkd=(5, 1, 3), label="5-qubit")
end

"""
    SteaneCode()

Construct Steane [7,1,3] code as a [`BasicCode`](@ref).
"""
function SteaneCode()
    BasicCode(["IIIXXXX", "IXXIIXX", "XIXIXIX", "IIIZZZZ", "IZZIIZZ", "ZIZIZIZ"],
        ["XXXXXXX"], ["ZZZZZZZ"]; nkd=(7, 1, 3), label="Steane")
end

end
