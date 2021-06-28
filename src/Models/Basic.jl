module Basic

import Qecsim.Model

using Qecsim.Model:StabilizerCode
using Qecsim:PauliTools as PT

struct BasicCode <: StabilizerCode
    stablizers::BitMatrix
    logical_xs::BitMatrix
    logical_zs::BitMatrix
    nkd::Tuple{Int, Int, Union{Int,Nothing}}
    label::String
end
function BasicCode(pauli_stabilizers, pauli_logical_xs, pauli_logical_zs;
                   nkd=nothing, label=nothing)
    nkd = !isnothing(nkd) ? nkd : (
        length(pauli_stabilizers) > 0 ? length(pauli_stabilizers[0]) : 0,
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

end
