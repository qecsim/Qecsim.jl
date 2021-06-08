"""
Tools for Pauli strings and binary symplectic vectors / matrices.
"""
module PauliTools

function hello(s)
    "Hello $s"
end

@doc raw"""
    bsp(A::BitArray, B::BitArray)

Return the binary symplectic product of A with B.

The binary symplectic product ``\odot`` is defined as
``A \odot B \equiv A \Lambda B \bmod 2`` where
``\Lambda = \left[\begin{smallmatrix} 0 & I \\ I & 0 \end{smallmatrix}\right]``.

# Examples
```jldoctest
using Qecsim: PauliTools as PT
a = BitVector([1, 0, 0, 0])  # XI
b = BitVector([0, 0, 1, 0])  # ZI
PT.bsp(a', b)
# output
1
```
```jldoctest
using Qecsim: PauliTools as PT
stabilizers = BitMatrix(  # 5-qubit stabilizers
    [1 0 0 1 0 0 1 1 0 0    # XZZXI
     0 1 0 0 1 0 0 1 1 0    # IXZZX
     1 0 1 0 0 0 0 0 1 1    # XIXZZ
     0 1 0 1 0 1 0 0 0 1])  # ZXIXZ
error = BitVector([0, 0, 1, 1, 0, 0, 1, 0, 1, 0])  # IZXYI
PT.bsp(stabilizers, error)
# output
4-element Vector{Int64}:
 0
 1
 1
 0
```
"""
function bsp(a::AbstractArray{Bool}, b::AbstractArray{Bool})
    # circshift b by half its 1st dimension to emulate symplectic product
    mod.(a * circshift(b, size(b, 1)/2), 2)  # mod elements to base 2
end

end