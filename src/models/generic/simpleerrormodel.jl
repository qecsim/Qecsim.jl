
"""
    SimpleErrorModel <: ErrorModel

Abstract supertype for simple IID error models that generate errors based on the number of
qubits and a probability distribution.
"""
abstract type SimpleErrorModel <: ErrorModel end

"""
    generate(error_model::SimpleErrorModel, code::StabilizerCode, p::Float64,
             [rng::AbstractRNG=GLOBAL_RNG]) -> BitVector

Generate a new IID error based on [`probability_distribution`](@ref). See also
[`Model.generate`](@ref).

!!! note

    The method [`probability_distribution`](@ref) should be implemented for concrete
    subtypes of [`SimpleErrorModel`](@ref).

"""
function Model.generate(error_model::SimpleErrorModel, code::StabilizerCode, p::Float64,
                        rng::AbstractRNG=GLOBAL_RNG)
    n_qubits = nkd(code)[1]
    weights = ProbabilityWeights(collect(probability_distribution(error_model, p)))
    error_pauli = String(sample(rng, ['I','X','Y','Z'], weights, n_qubits))
    return to_bsf(error_pauli)
end

struct DepolarizingErrorModel <: SimpleErrorModel end
function Model.label(::DepolarizingErrorModel)
    return "Depolarizing"
end
function Model.probability_distribution(::DepolarizingErrorModel, p::Float64)
    px = py = pz = p / 3
    pi = 1 - sum((px, py, pz))
    return (pi, px, py, pz)
end

struct BitFlipErrorModel <: SimpleErrorModel end
function Model.label(::BitFlipErrorModel)
    return "Bit-flip"
end
function Model.probability_distribution(::BitFlipErrorModel, p::Float64)
    return (1-p, p, 0, 0)  # (pi, px, py, pz)
end

struct PhaseFlipErrorModel <: SimpleErrorModel end
function Model.label(::PhaseFlipErrorModel)
    return "Phase-flip"
end
function Model.probability_distribution(::PhaseFlipErrorModel, p::Float64)
    return (1-p, 0, 0, p)  # (pi, px, py, pz)
end

struct BitPhaseFlipErrorModel <: SimpleErrorModel end
function Model.label(::BitPhaseFlipErrorModel)
    return "Bit-phase-flip"
end
function Model.probability_distribution(::BitPhaseFlipErrorModel, p::Float64)
    return (1-p, 0, p, 0)  # (pi, px, py, pz)
end
