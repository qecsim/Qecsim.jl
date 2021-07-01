"""
Tools for Pauli strings and binary symplectic vectors / matrices.
"""
module PauliTools

export bsp, to_pauli, to_bsf

@doc raw"""
    bsp(A::Union{BitVector, BitMatrix}, B::Union{BitVector, BitMatrix})

Return the binary symplectic product of A with B, given in binary symplectic form.

The binary symplectic product ``\odot`` is defined as
``A \odot B \equiv A \Lambda B \bmod 2`` where
``\Lambda = \left[\begin{smallmatrix} 0 & I \\ I & 0 \end{smallmatrix}\right]``.

# Examples
```jldoctest
julia> a = BitVector([1, 0, 0, 0]);  # XI

julia> b = BitVector([0, 0, 1, 0]);  # ZI

julia> bsp(a', b)
1
```
```jldoctest
julia> stabilizers = BitMatrix(  # 5-qubit stabilizers
       [1 0 0 1 0 0 1 1 0 0     # XZZXI
        0 1 0 0 1 0 0 1 1 0     # IXZZX
        1 0 1 0 0 0 0 0 1 1     # XIXZZ
        0 1 0 1 0 1 0 0 0 1]);  # ZXIXZ

julia> error = BitVector([0, 0, 1, 1, 0, 0, 1, 0, 1, 0]);  # IZXYI

julia> bsp(stabilizers, error)
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
    to_pauli(bsf::Union{BitVector, BitMatrix})

Convert binary symplectic form to Pauli string operator(s).

A vector is converted to a single Pauli string. A matrix is converted row-by-row to a
collection of Pauli strings.

# Examples
```jldoctest
julia> to_pauli(BitVector([1, 0, 0, 0, 1, 0, 0, 1, 0, 1]))
"XIZIY"
```
```jldoctest
julia> to_pauli(BitMatrix([1 0 0 0 1 0 0 1 0 1; 0 1 0 1 0 0 0 1 1 0]))
2-element Vector{String}:
 "XIZIY"
 "IXZYI"
```
"""
function to_pauli(bsf::AbstractVector{Bool})
    l = length(bsf)  # bsf = (1 0 0 0 1 | 0 0 1 0 1)
    x, z = bsf[1:l÷2], bsf[l÷2+1:end]  # x = (1 0 0 0 1), z = (0 0 1 0 1)
    String(map(i -> "IXZY"[i+1], x+2z))  # x+2z = (1 0 2 0 3) -> "XIZIY"
end
function to_pauli(bsfs::AbstractMatrix{Bool})
    String[to_pauli(bsf) for bsf in eachrow(bsfs)]
end

"""
    to_bsf(pauli::Union{String, Iterable of String})

Convert the Pauli string operator(s) to binary symplectic form.

A single Pauli string is converted to a vector. A collection of Pauli strings is converted
to a matrix where each row corresponds to a Pauli.

# Examples
```jldoctest
julia> to_bsf("XIZIY")
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
julia> to_bsf(["XIZIY", "IXZYI"])
2×10 BitMatrix:
 1  0  0  0  1  0  0  1  0  1
 0  1  0  1  0  0  0  1  1  0
```
"""
function to_bsf(pauli::AbstractString)
    p = collect(pauli)
    ys = p .== 'Y'
    vcat((p .== 'X') .| ys, (p .== 'Z') .| ys)
end
function to_bsf(paulis)
    # Note: I cannot see good way to restrict this method based on eltype of AbstractString
    # so an 'isa' type test is included to fail-fast and avoid infinite recursion.

    # for each pauli, if string convert to bsf vector, transpose and concatenate as rows
    vcat((isa(p, AbstractString) ? transpose(to_bsf(p))
          : throw(ArgumentError("invalid Pauli type: $(typeof(p))"))
          for p in paulis)...)
end

end
