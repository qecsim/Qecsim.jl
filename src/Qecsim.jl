"""
Package for simulating quantum error correction using stabilizer codes.
"""
module Qecsim

# imports
using Reexport

# exports
export QecsimError
include("error.jl")

# include core sub-modules (and reexport)
include("PauliTools.jl")
include("Model.jl")
include("App.jl")
@reexport using .PauliTools
@reexport using .Model
@reexport using .App

# include implementation sub-modules (not reexported)
include("models/BasicModels.jl")
include("models/GenericModels.jl")

end
