
"""
    SimpleErrorModel <: ErrorModel

Abstract supertype for simple IID error models that generate errors based on the number of
qubits and a probability distribution.
"""
abstract type SimpleErrorModel <: ErrorModel end

"""
    generate(error_model::SimpleErrorModel, code::StabilizerCode, p::Float64,
             [rng::AbstractRNG=GLOBAL_RNG]) -> BitVector

Generate a new IID error based on [`Model.probability_distribution`](@ref). See also
[`Model.generate`](@ref).

!!! note

    The method [`Model.probability_distribution`](@ref) should be implemented for concrete
    subtypes of [`SimpleErrorModel`](@ref).

"""
function Model.generate(error_model::SimpleErrorModel, code::StabilizerCode, p::Float64,
                        rng::AbstractRNG=GLOBAL_RNG)
    n_qubits = nkd(code)[1]
    weights = ProbabilityWeights(collect(probability_distribution(error_model, p)))
    error_pauli = String(sample(rng, ['I','X','Y','Z'], weights, n_qubits))
    return to_bsf(error_pauli)
end

@doc raw"""
    BitFlipErrorModel <: SimpleErrorModel

IID error model with probability vector: ``(p_I, p_X, p_Y p_Z) = (1-p, p, 0, 0)``,
where ``p`` is the probability of an error on a single-qubit.
"""
struct BitFlipErrorModel <: SimpleErrorModel end
Model.label(::BitFlipErrorModel) = "Bit-flip"
Model.probability_distribution(::BitFlipErrorModel, p::Float64) = (1 - p, p, 0., 0.)

@doc raw"""
    BitPhaseFlipErrorModel <: SimpleErrorModel

IID error model with probability vector: ``(p_I, p_X, p_Y p_Z) = (1-p, 0, p, 0)``,
where ``p`` is the probability of an error on a single-qubit.
"""
struct BitPhaseFlipErrorModel <: SimpleErrorModel end
Model.label(::BitPhaseFlipErrorModel) = "Bit-phase-flip"
Model.probability_distribution(::BitPhaseFlipErrorModel, p::Float64) = (1 - p, 0., p, 0.)

@doc raw"""
    DepolarizingErrorModel <: SimpleErrorModel

IID error model with probability vector: ``(p_I, p_X, p_Y p_Z) = (1-p, p/3, p/3, p/3)``,
where ``p`` is the probability of an error on a single-qubit.
"""
struct DepolarizingErrorModel <: SimpleErrorModel end
Model.label(::DepolarizingErrorModel) = "Depolarizing"
function Model.probability_distribution(::DepolarizingErrorModel, p::Float64)
    px = py = pz = p / 3
    pi = 1 - sum((px, py, pz))
    return (pi, px, py, pz)
end

@doc raw"""
    PhaseFlipErrorModel <: SimpleErrorModel

IID error model with probability vector: ``(p_I, p_X, p_Y p_Z) = (1-p, 0, 0, p)``,
where ``p`` is the probability of an error on a single-qubit.
"""
struct PhaseFlipErrorModel <: SimpleErrorModel end
Model.label(::PhaseFlipErrorModel) = "Phase-flip"
Model.probability_distribution(::PhaseFlipErrorModel, p::Float64) = (1 - p, 0., 0., p)
