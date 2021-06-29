"""
Tools for Pauli strings and binary symplectic vectors / matrices.
"""
module PauliTools

@doc raw"""
    bsp(A::Union{BitVector, BitMatrix}, B::Union{BitVector, BitMatrix})

Return the binary symplectic product of A with B, given in binary symplectic form.

The binary symplectic product ``\odot`` is defined as
``A \odot B \equiv A \Lambda B \bmod 2`` where
``\Lambda = \left[\begin{smallmatrix} 0 & I \\ I & 0 \end{smallmatrix}\right]``.

# Examples
```jldoctest
julia> using Qecsim: PauliTools as PT

julia> a = BitVector([1, 0, 0, 0]);  # XI

julia> b = BitVector([0, 0, 1, 0]);  # ZI

julia> PT.bsp(a', b)
1
```
```jldoctest
julia> using Qecsim: PauliTools as PT

julia> stabilizers = BitMatrix(  # 5-qubit stabilizers
       [1 0 0 1 0 0 1 1 0 0     # XZZXI
        0 1 0 0 1 0 0 1 1 0     # IXZZX
        1 0 1 0 0 0 0 0 1 1     # XIXZZ
        0 1 0 1 0 1 0 0 0 1]);  # ZXIXZ

julia> error = BitVector([0, 0, 1, 1, 0, 0, 1, 0, 1, 0]);  # IZXYI

julia> PT.bsp(stabilizers, error)
4-element Vector{Int64}:
 0
 1
 1
 0
```
"""
function bsp(a::AbstractVecOrMat{Bool}, b::AbstractVecOrMat{Bool})
    # circshift b by half its 1st dimension to emulate symplectic product
    mod.(a * circshift(b, size(b, 1)/2), 2)  # mod elements to base 2
end

"""
    pauli_to_bsf(pauli)

Convert the Pauli operator(s) to binary symplectic form.

A single string is converted to a vector. A collection of strings is converted to a matrix
where each row corresponds to a pauli.

# Examples
```jldoctest
julia> using Qecsim: PauliTools as PT

julia> PT.pauli_to_bsf("XIZIY")
10-element BitVector:
 1
 0
 0
 0
 1
 0
 0
 1
 0
 1
```
```jldoctest
julia> using Qecsim: PauliTools as PT

julia> PT.pauli_to_bsf(["XIZIY", "IXZYI"])
2×10 BitMatrix:
 1  0  0  0  1  0  0  1  0  1
 0  1  0  1  0  0  0  1  1  0
```
"""
function pauli_to_bsf(pauli::AbstractString)
    p = collect(pauli)
    ys = p .== 'Y'
    vcat((p .== 'X') .| ys, (p .== 'Z') .| ys)
end
function pauli_to_bsf(paulis)
    vcat((transpose(pauli_to_bsf(p)) for p in paulis)...)
end

function bsf_to_pauli(bsf::AbstractVector{Bool})
    nothing
end
function bsf_to_pauli(bsfs::AbstractMatrix{Bool})
    nothing
end

end
