"""
Package for simulating quantum error correction using stabilizer codes.
"""
module Qecsim

using Reexport

export doubler
include("poc.jl")

# include core sub-modules
include("PauliTools.jl")
include("Model.jl")
include("Models/Models.jl")

# reexport core sub-module exports
@reexport using Qecsim.PauliTools
@reexport using Qecsim.Model

end
