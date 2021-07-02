"""
Abstract types and methods for codes, error models and decoders.
"""
module Model

using LinearAlgebra:I
using Qecsim:QecsimError
using Qecsim.PauliTools:bsp

export StabilizerCode
export stabilizers, logical_xs, logical_zs, logicals, nkd, label, validate

"""
Abstract supertype for stabilizer codes.
"""
abstract type StabilizerCode end

"""
    stabilizers(code::StabilizerCode) -> BitMatrix

Return the stabilizers in binary symplectic form.

Each row is a stabilizer generator. An overcomplete set of generators can be included to
simplify decoding.

!!! note "Abstract method"

    This method should be implemented for concrete subtypes.
"""
function stabilizers end

"""
    logical_xs(code::StabilizerCode) -> BitMatrix

Return the logical X operators in binary symplectic form.

Each row is an operator. The order should match that of [`logical_zs`](@ref).

!!! note "Abstract method"

    This method should be implemented for concrete subtypes.
"""
function logical_xs end

"""
    logical_zs(code::StabilizerCode) -> BitMatrix

Return the logical Z operators in binary symplectic form.

Each row is an operator. The order should match that of [`logical_xs`](@ref).

!!! note "Abstract method"

    This method should be implemented for concrete subtypes.
"""
function logical_zs end

"""
    logicals(code::StabilizerCode) -> BitMatrix

Return the logical operators in binary symplectic form.

Each row is an operator. X operators are stacked above Z operators in the order given by
[`logical_xs`](@ref) and [`logical_zs`](@ref).
"""
function logicals(code::StabilizerCode)
    vcat(logical_xs(code), logical_zs(code))
end

"""
    nkd(code::StabilizerCode) -> Tuple{Int, Int, Union{Int,Nothing}}

Return a descriptor in the format `(n, k, d)`, where `n` is the number of physical qubits,
`k` is the number of logical qubits, and `d` is the distance of the code (or `nothing` if
unknown).

!!! note "Abstract method"

    This method should be implemented for concrete subtypes.
"""
function nkd end

"""
    label(code::StabilizerCode) -> String

Return a label suitable for use in plots and for grouping results.

!!! note "Abstract method"

    This method should be implemented for concrete subtypes.
"""
function label end

@doc raw"""
    validate(code::StabilizerCode)

Perform sanity checks.

If any of the following fail then a [`QecsimError`](@ref) is thrown:

  * ``S \odot S^T = 0``
  * ``S \odot L^T = 0``
  * ``L \odot L^T = \Lambda``

where ``S`` and ``L`` are the code [`stabilizers`](@ref) and [`logicals`](@ref),
respectively, and ``\odot`` and ``\Lambda`` are defined in [`bsp`](@ref).
"""
function validate(code::StabilizerCode)
    s, l = stabilizers(code), logicals(code)
    all(bsp(s, transpose(s)) .== 0) || throw(QecsimError(
        "stabilizers do not mutually commute"))
    all(bsp(s, transpose(l)) .== 0) || throw(QecsimError(
        "stabilizers do not commute with logicals"))
    # twisted identity with same size as logicals
    nlogicals = size(l, 1)
    twistedI = circshift(Matrix(I, nlogicals, nlogicals), nlogicals / 2)
    (bsp(l, transpose(l)) == twistedI) || throw(QecsimError(
        "logicals do not mutually twist commute"))
    return
end

end
