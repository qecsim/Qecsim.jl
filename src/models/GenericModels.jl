"""
Generic error models and decoders compatible with any stabilizer codes.
"""
module GenericModels

export SimpleErrorModel
export BitFlipErrorModel
export BitPhaseFlipErrorModel
export DepolarizingErrorModel
export PhaseFlipErrorModel

import Qecsim.Model

using Qecsim.Model
using Qecsim.PauliTools:to_bsf
using Random:AbstractRNG, GLOBAL_RNG
using StatsBase:ProbabilityWeights, sample

include("generic/simpleerrormodel.jl")

end
