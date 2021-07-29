"""
Basic stabilizer codes
"""
module BasicModels

# imports
using ..Model
using ..PauliTools:to_bsf, to_pauli

# exports
export BasicCode, FiveQubitCode, SteaneCode

struct BasicCode <: StabilizerCode
    stabilizers::BitMatrix
    logical_xs::BitMatrix
    logical_zs::BitMatrix
    nkd::Tuple{Int,Int,Union{Int,Missing}}
    label::String
end
"""
    BasicCode <: StabilizerCode

    BasicCode(pauli_stabilizers::AbstractVector{<:AbstractString},
              pauli_logical_xs::AbstractVector{<:AbstractString},
              pauli_logical_zs::AbstractVector{<:AbstractString},
              nkd::Tuple{Int,Int,Union{Int,Missing}}=nothing, label::AbstractString=nothing)

Construct a basic code from string representations of stabilizers and logical operators.

Paulis are expressed as strings of capitalized I, X, Y, Z characters, with one character per
physical qubit. Logical X and Z operators are in matching order, with one of each for each
logical qubit. Optional `nkd` defaults to `n` and `k` evaluated and `d` missing. Optional
`label` defaults to "Basic [n,k,d]".

# Examples
```jldoctest
julia> using Qecsim.BasicModels

julia> code = BasicCode(["ZZI", "IZZ"], ["XXX"], ["IIZ"]);  # 3-qubit repetition

julia> validate(code)  # no error indicates operators satisfy commutation relations

julia> nkd(code)  # default nkd
(3, 1, missing)

julia> label(code)  # default label
"Basic [3,1,missing]"
```
"""
function BasicCode(pauli_stabilizers::AbstractVector{<:AbstractString},
                    pauli_logical_xs::AbstractVector{<:AbstractString},
                    pauli_logical_zs::AbstractVector{<:AbstractString},
                    nkd=nothing, label=nothing)
    nkd = !isnothing(nkd) ? nkd : (
        length(pauli_stabilizers) > 0 ? length(pauli_stabilizers[1]) : 0,
        length(pauli_logical_xs),
        missing
    )
    label = !isnothing(label) ? label : "Basic [$(nkd[1]),$(nkd[2]),$(nkd[3])]"
    return BasicCode(to_bsf(pauli_stabilizers), to_bsf(pauli_logical_xs),
                     to_bsf(pauli_logical_zs), nkd, label)
end
Model.stabilizers(code::BasicCode) = code.stabilizers
Model.logical_xs(code::BasicCode) = code.logical_xs
Model.logical_zs(code::BasicCode) = code.logical_zs
Model.nkd(code::BasicCode) = code.nkd
Model.label(code::BasicCode) = code.label
function Base.show(io::IO, x::BasicCode)
    show(io, typeof(x)); print(io, '(')
    for f in (x.stabilizers, x.logical_xs, x.logical_zs)
        show(io, to_pauli(f)); print(io, ", ")
    end
    show(io, x.nkd); print(io, ", ")
    show(io, x.label); print(io, ')')
end

"""
    FiveQubitCode() -> BasicCode

Construct 5-qubit [5,1,3] code as a [`BasicCode`](@ref).
"""
function FiveQubitCode()
    return BasicCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"], ["XXXXX"], ["ZZZZZ"], (5, 1, 3),
                      "5-qubit")
end

"""
    SteaneCode() -> BasicCode

Construct Steane [7,1,3] code as a [`BasicCode`](@ref).
"""
function SteaneCode()
    return BasicCode(["IIIXXXX", "IXXIIXX", "XIXIXIX", "IIIZZZZ", "IZZIIZZ", "ZIZIZIZ"],
                     ["XXXXXXX"], ["ZZZZZZZ"], (7, 1, 3), "Steane")
end

end
