"""
    NaiveDecoder <: Decoder

Naive decoder that iterates through all possible errors, in ascending weight, and resolves
to the first error that matches the syndrome.

!!! note

    This decoder is slow for large number of qubits and high weight errors.
"""
struct NaiveDecoder <: Decoder end
Model.label(::NaiveDecoder) = "Naive"
function Model.decode(::NaiveDecoder, code::StabilizerCode, syndrome::AbstractVector{Bool};
                      kwargs...)
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
