"""
Tools for Pauli strings and binary symplectic vectors / matrices.
"""
module PauliTools

# exports
export bsp, pack, unpack, to_bsf, to_pauli, weight

@doc raw"""
    bsp(A::AbstractVecOrMat{Bool}, B::AbstractVecOrMat{Bool})
        -> Union{Bool,BitVector,BitMatrix}

Return the binary symplectic product of `A` with `B`, given in binary symplectic form.

The binary symplectic product ``\odot`` is defined as
``A \odot B \equiv A \Lambda B \bmod 2`` where
``\Lambda = \left[\begin{smallmatrix} 0 & I \\ I & 0 \end{smallmatrix}\right]``.

# Examples
```jldoctest
julia> a = BitVector([1, 0, 0, 0]);  # XI

julia> b = BitVector([0, 0, 1, 0]);  # ZI

julia> bsp(a', b)
true
```
```jldoctest
julia> stabilizers = BitMatrix(  # 5-qubit stabilizers
       [1 0 0 1 0 0 1 1 0 0     # XZZXI
        0 1 0 0 1 0 0 1 1 0     # IXZZX
        1 0 1 0 0 0 0 0 1 1     # XIXZZ
        0 1 0 1 0 1 0 0 0 1]);  # ZXIXZ

julia> error = BitVector([0, 0, 1, 1, 0, 0, 1, 0, 1, 0]);  # IZXYI

julia> bsp(stabilizers, error)
4-element BitVector:
 0
 1
 1
 0
```
"""
function bsp(a::AbstractVecOrMat{Bool}, b::AbstractVecOrMat{Bool})
    # circshift b by half its 1st dimension to emulate symplectic product
    return isodd.(a * circshift(b, size(b, 1)/2))  # mod elements to base 2
end

"""
    pack(bsf::AbstractVector{Bool}) -> Tuple{String,Int}

Pack a binary vector into a concise representation, typically for log output. See also
[`unpack`](@ref).

# Examples
```jldoctest
julia> a = BitVector([1, 0, 1, 0, 1, 1]);  # XZY

julia> b = pack(a)  # (hex_value, length)
("2b", 6)
julia> unpack(b) == a
true
```
"""
function pack(bsf::AbstractVector{Bool})
    bit_str = join(d ? '1' : '0' for d in bsf)  # convert to bit-string
    val = parse(BigInt, bit_str; base=2)  # parse to val (BigInt to support long vectors)
    return string(val, base=16), length(bsf)  # hex_str, len
end

"""
    unpack(packed_bsf::Tuple{String,Int}) -> BitVector

Unpack a binary vector from a concise representation, typically from log output. See also
[`pack`](@ref).

# Examples
```jldoctest
julia> a = ("2b", 6);  # (hex_value, length)

julia> b = unpack(a)  # XZY
6-element BitVector:
 1
 0
 1
 0
 1
 1
julia> pack(b) == a
true
```
"""
function unpack(packed_bsf::Tuple{String,Int})
    hex_str, len = packed_bsf
    val = parse(BigInt, hex_str; base=16)  # parse to val (BigInt to support long vectors)
    return reverse!(BitVector(digits(val; base=2, pad=len)))  # binary-digits as BitVector
end

"""
    to_bsf(pauli::Union{AbstractString,AbstractVector{<:AbstractString}})
        -> Union{BitVector,BitMatrix}

Convert the Pauli string operator(s) to binary symplectic form.

A single Pauli string is converted to a vector. A vector of Pauli strings is converted to a
matrix where each row corresponds to a Pauli.

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
    return vcat((p .== 'X') .| ys, (p .== 'Z') .| ys)
end
function to_bsf(paulis::AbstractVector{<:AbstractString})
    # for each pauli string, convert to bsf vector, transpose; then concatenate as rows
    return vcat((transpose(to_bsf(p)) for p in paulis)...)
end

"""
    to_pauli(bsf::AbstractVecOrMat{Bool}) -> Union{String,Vector{String}}

Convert the binary symplectic form to Pauli string operator(s).

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
    return String(map(i -> "IXZY"[i+1], x+2z))  # x+2z = (1 0 2 0 3) -> "XIZIY"
end
function to_pauli(bsfs::AbstractMatrix{Bool})
    return String[to_pauli(bsf) for bsf in eachrow(bsfs)]
end

"""
    weight(bsf::AbstractVecOrMat{Bool}) -> Int

Return the weight of the binary symplectic form.

# Examples
```jldoctest
julia> weight(BitVector([1, 0, 0, 0, 1, 0, 0, 1, 0, 1]))
3
```
```jldoctest
julia> weight(BitMatrix([1 0 0 0 1 0 0 1 0 1; 1 1 1 1 1 0 0 0 0 0]))
8
```
"""
function weight(bsf::AbstractVector{Bool})
    l = length(bsf)  # bsf = (1 0 0 | 1 1 0)
    return count(>(0), bsf[1:l÷2] + bsf[l÷2+1:end])  # count(>(0), (2 1 0))
end
function weight(bsfs::AbstractMatrix{Bool})
    l = size(bsfs, 2)  # bsf = (1 0 0 | 1 1 0; 1 1 1 | 0 1 0)
    return count(>(0), bsfs[:,1:l÷2] + bsfs[:,l÷2+1:end])  # count(>(0), (2 1 0; 1 2 1))
end

"""
    weight(pauli::Union{AbstractString,AbstractVector{<:AbstractString}}) -> Int

Return the weight of the Pauli string operator(s).

# Examples
```jldoctest
julia> weight("XIZIY")
3
```
```jldoctest
julia> weight(["XIZIY", "XXXXX"])
8
```
"""
function weight(pauli::AbstractString)
    return count(!=('I'), pauli)  # count(!=('I'), "IXYIZ")
end
function weight(paulis::AbstractVector{<:AbstractString})
    return sum(weight, paulis)
end

end
