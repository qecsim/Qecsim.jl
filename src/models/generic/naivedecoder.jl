"""
    NaiveDecoder <: Decoder

Naive decoder that iterates through all possible errors, in ascending weight, and resolves
to the first error that matches the syndrome.

!!! note

    This decoder is slow for even moderate numbers of qubits. By default, it is restricted
    to codes with a maximum of 10 qubits.
"""
struct NaiveDecoder <: Decoder
    max_qubits::Int
    NaiveDecoder(max_qubits=10) = new(max_qubits)
end
Model.label(::NaiveDecoder) = "Naive"
function Model.decode(
    decoder::NaiveDecoder,
    code::StabilizerCode,
    syndrome::AbstractVector{Bool};
    kwargs...
)
    if 0 <= decoder.max_qubits < nkd(code)[1]
        throw(ArgumentError("NaiveDecoder restricted to $(decoder.max_qubits) qubits"))
    end
    recovery = _minimum_weight_recovery(code, syndrome)
    return DecodeResult(recovery)
end

function _minimum_weight_recovery(code, syndrome)
    n_qubits = nkd(code)[1]
    min_weight = 0
    max_weight = n_qubits
    # iterate recovery operations of increasing weight
    for weight in min_weight:max_weight
        for qubit_indices in combinations(1:n_qubits, weight)
            for xyzs in Iterators.product(fill(1:3, weight)...)
                # pauli encoded as 0123 for IXZY
                v = zeros(Int, n_qubits)
                v[qubit_indices] = collect(xyzs)
                # pauli encoded in binary symplectic form
                recovery = BitVector(vcat(v .% 2, v .÷ 2))
                if bsp(stabilizers(code), recovery) == syndrome
                    return recovery
                end
            end
        end
    end
    return falses(n_qubits * 2)
end
