"""
Package for simulating quantum error correction using stabilizer codes.
"""
module Qecsim

export doubler
include("poc.jl")

include("PauliTools.jl")
include("Model.jl")
include("Models/Models.jl")

end
