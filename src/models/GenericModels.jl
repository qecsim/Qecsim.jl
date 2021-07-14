"""
Generic error models and decoders compatible with any stabilizer codes.
"""
module GenericModels

# exports
#  error models
export SimpleErrorModel
export BitFlipErrorModel
export BitPhaseFlipErrorModel
export DepolarizingErrorModel
export PhaseFlipErrorModel
#  decoders
export NaiveDecoder

# imports
import Qecsim.Model
# usings
using Qecsim.Model
using Qecsim.PauliTools: bsp, to_bsf
using Combinatorics: combinations
using Random: AbstractRNG, GLOBAL_RNG
using StatsBase: ProbabilityWeights, sample

# includes
include("generic/simpleerrormodel.jl")
include("generic/naivedecoder.jl")

end
