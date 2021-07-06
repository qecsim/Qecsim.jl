"""
Package for simulating quantum error correction using stabilizer codes.
"""
module Qecsim

using Reexport

export QecsimError
include("error.jl")

# include core sub-modules (and reexport)
include("PauliTools.jl")
include("Model.jl")
@reexport using Qecsim.PauliTools
@reexport using Qecsim.Model

# include implementation sub-modules (not reexported)
include("models/BasicModels.jl")
include("models/GenericModels.jl")

end
