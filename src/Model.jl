"""
Abstract types and methods for codes, error models and decoders.
"""
module Model

# imports
using ..Qecsim:QecsimError
using ..PauliTools:bsp
using LinearAlgebra:I

# exports
# - abstract model
export AbstractModel
export label
# - stabilizer code
export StabilizerCode
export stabilizers, logical_xs, logical_zs, logicals, nkd, validate
# - error model
export ErrorModel
export generate, probability_distribution
# - decoder
export Decoder
export DecodeResult
export decode


# AbstractModel

"""
Abstract supertype for models.
"""
abstract type AbstractModel end

"""
    label(model) -> String

Return a label suitable for use in plots and for grouping results.

!!! note "Abstract method"

    This method should be implemented for concrete subtypes or duck-typed implementations of
    [`AbstractModel`](@ref).
"""
function label end


# StabilizerCode

"""
    StabilizerCode <: AbstractModel

Abstract supertype for stabilizer codes.
"""
abstract type StabilizerCode <: AbstractModel end

"""
    stabilizers(code) -> BitMatrix

Return the stabilizers in binary symplectic form.

Each row is a stabilizer generator. An overcomplete set of generators can be included to
simplify decoding.

!!! note "Abstract method"

    This method should be implemented for concrete subtypes or duck-typed implementations of
    [`StabilizerCode`](@ref).
"""
function stabilizers end

"""
    logical_xs(code) -> BitMatrix

Return the logical X operators in binary symplectic form.

Each row is an operator. The order should match that of [`logical_zs`](@ref).

!!! note "Abstract method"

    This method should be implemented for concrete subtypes or duck-typed implementations of
    [`StabilizerCode`](@ref).
"""
function logical_xs end

"""
    logical_zs(code) -> BitMatrix

Return the logical Z operators in binary symplectic form.

Each row is an operator. The order should match that of [`logical_xs`](@ref).

!!! note "Abstract method"

    This method should be implemented for concrete subtypes or duck-typed implementations of
    [`StabilizerCode`](@ref).
"""
function logical_zs end

"""
    logicals(code) -> BitMatrix

Return the logical operators in binary symplectic form.

Each row is an operator. X operators are stacked above Z operators in the order given by
[`logical_xs`](@ref) and [`logical_zs`](@ref).
"""
function logicals(code)
    return vcat(logical_xs(code), logical_zs(code))
end

"""
    nkd(code) -> Tuple{Int,Int,Union{Int,Missing}}

Return a descriptor in the format `(n, k, d)`, where `n` is the number of physical qubits,
`k` is the number of logical qubits, and `d` is the distance of the code (or `missing` if
unknown).

!!! note "Abstract method"

    This method should be implemented for concrete subtypes or duck-typed implementations of
    [`StabilizerCode`](@ref).
"""
function nkd end

@doc raw"""
    validate(code)

Perform sanity checks.

If any of the following fail then a [`QecsimError`](@ref) is thrown:

* ``S ⊙ S^T = 0``
* ``S ⊙ L^T = 0``
* ``L ⊙ L^T = Λ``

where ``S`` and ``L`` are the code [`stabilizers`](@ref) and [`logicals`](@ref),
respectively, and ``⊙`` and ``Λ`` are defined in [`PauliTools.bsp`](@ref bsp).
"""
function validate(code)
    s, l = stabilizers(code), logicals(code)
    !any(bsp(s, transpose(s))) || throw(QecsimError("stabilizers do not mutually commute"))
    !any(bsp(s, transpose(l))) || throw(QecsimError(
        "stabilizers do not commute with logicals"))
    # twisted identity with same size as logicals
    nlogicals = size(l, 1)
    twistedI = circshift(Matrix(I, nlogicals, nlogicals), nlogicals / 2)
    (bsp(l, transpose(l)) == twistedI) || throw(QecsimError(
        "logicals do not mutually twist commute"))
    return nothing
end


# ErrorModel

"""
    ErrorModel <: AbstractModel

Abstract supertype for error models.
"""
abstract type ErrorModel <: AbstractModel end

"""
    generate(error_model, code, p::Real, [rng::AbstractRNG=GLOBAL_RNG]) -> BitVector

Generate a new error in binary symplectic form according to the `error_model` and `code`,
where `p` is typically the probability of an error on a single qubit.

!!! note "Abstract method"

    This method should be implemented for concrete subtypes or duck-typed implementations of
    [`ErrorModel`](@ref).
"""
function generate end

"""
    probability_distribution(error_model, p::Real) -> NTuple{4,Real}

Return the single-qubit probability distribution amongst Pauli I, X, Y and Z, where `p` is
the overall probability of an error on a single qubit.

!!! note "Abstract method [optional]"

    This method is **not** invoked by any core modules. Since it is often useful for
    decoders, it is provided as a template and concrete subtypes or duck-typed
    implementations of [`ErrorModel`](@ref) are encouraged to implement it when appropriate,
    particularly for IID error models.
"""
function probability_distribution end


# Decoder
"""
    Decoder <: AbstractModel

Abstract supertype for decoders.
"""
abstract type Decoder <: AbstractModel end

"""
    decode(decoder, code, syndrome::AbstractVector{Bool}; kwargs...) -> DecodeResult

Resolve a recovery operation for the given `code` and `syndrome`, or evaluate the success of
decoding, as encapsulated in the decode result.

The syndrome has length equal to the number of code stabilizers, and element values of 0
or 1 (false or true) indicate whether the corresponding stabilizer does or does not commute
with the error, respectively.

Keyword parameters `kwargs` may be provided by the client, e.g. [`App`](@ref App), with
context values such as `error_model`, `error_probability` and `error`. Most implementations
will ignore such parameters; however, if they are used, implementations should declare them
explicitly and treat them as optional.

See also [`DecodeResult`](@ref).

!!! note "Abstract method"

    This method should be implemented for concrete subtypes or duck-typed implementations of
    [`Decoder`](@ref).
"""
function decode end

"""
    DecodeResult(
        success::Union{Nothing,Bool},
        recovery::Union{Nothing,AbstractVector{Bool}},
        logical_commutations::Union{Nothing,AbstractVector{Bool}}
        custom_values::Union{Nothing,AbstractVector}
    )
    DecodeResult(;
        success::Union{Nothing,Bool}=nothing,
        recovery::Union{Nothing,AbstractVector{Bool}}=nothing,
        logical_commutations::Union{Nothing,AbstractVector{Bool}}=nothing
        custom_values::Union{Nothing,AbstractVector}=nothing
    )

Construct a decoding result as returned by [`decode`](@ref).

Typically decoders will provide a `recovery` operation and delegate the evaluation of
`success` and `logical_commutations` to the client, e.g. [`App`](@ref App). Optionally,
`success` and/or `logical_commutations` may be provided as overrides. At least one of
`recovery` or `success` must be specified to allow a success value to be resolved.
Additionally `custom_values` may be specified. If `logical_commutations` and/or
`custom_values` are provided then they should be of consistent type and size over
identically parameterized simulation runs. For example, [`App`](@ref App) will sum
`logical_commutations` across runs and, similarly, if `custom_values` are numbers they will
be summed across runs, and if they are arrays they will be concatenated using `vcat`.

See also [`decode`](@ref).

# Examples
```jldoctest
julia> DecodeResult(; recovery=BitVector([1, 0, 1, 0, 1, 1]))  # typical use-case
DecodeResult{Nothing}(nothing, Bool[1, 0, 1, 0, 1, 1], nothing, nothing)

julia> DecodeResult(; success=true)  # override success
DecodeResult{Nothing}(true, nothing, nothing, nothing)

julia> DecodeResult(; success=false, logical_commutations=BitVector([1, 0]))  # override all
DecodeResult{Nothing}(false, nothing, Bool[1, 0], nothing)

julia> DecodeResult(; success=true, custom_values=[2.3, 4.1])  # custom values (numbers)
DecodeResult{Vector{Float64}}(true, nothing, nothing, [2.3, 4.1])

julia> DecodeResult(; success=true, custom_values=[[2.3], [4.1]])  # custom values (arrays)
DecodeResult{Vector{Vector{Float64}}}(true, nothing, nothing, [[2.3], [4.1]])

julia> DecodeResult(true, nothing, nothing, [2.3, 4.1])  # positional parameters
DecodeResult{Vector{Float64}}(true, nothing, nothing, [2.3, 4.1])

julia> DecodeResult(; success=nothing, recovery=nothing)  # too few specified parameters
ERROR: QecsimError: at least one of 'success' or 'recovery' must be specified
Stacktrace:
[...]
```
"""
struct DecodeResult{T<:Union{Nothing,Vector}}
    # Parametric type to avoid abstract types in struct, see:
    # https://docs.julialang.org/en/v1/manual/performance-tips/#Type-declarations
    # https://discourse.julialang.org/t/union-of-nothing-and-parametric-types/11986/5
    success::Union{Nothing,Bool}
    recovery::Union{Nothing,BitVector}
    logical_commutations::Union{Nothing,BitVector}
    custom_values::T
    # lenient constructor that infers types
    function DecodeResult(success, recovery, logical_commutations,
                          custom_values::Union{Nothing,AbstractVector{S}}) where S
        if isnothing(success) && isnothing(recovery)
            throw(QecsimError("at least one of 'success' or 'recovery' must be specified"))
        end
        if isnothing(custom_values)
            new{Nothing}(success, recovery, logical_commutations, custom_values)
        else
            new{Vector{S}}(success, recovery, logical_commutations, collect(custom_values))
        end
    end
end
# keyword convenience constructor
function DecodeResult(; success=nothing, recovery=nothing, logical_commutations=nothing,
                      custom_values=nothing)
    DecodeResult(success, recovery, logical_commutations, custom_values)
end

end
