using Test
using SafeTestsets

@testset verbose=true "Qecsim.jl" begin
    @safetestset "error.jl" begin include("error.jl") end
    @safetestset "PauliTools.jl" begin include("PauliTools.jl") end
    @safetestset "Model.jl" begin include("Model.jl") end
    @safetestset "models/BasicModels.jl" begin include("models/BasicModels.jl") end
    @safetestset "models/GenericModels/simpleerrormodel.jl" begin
        include("models/generic/simpleerrormodel.jl")
    end
    @safetestset "models/GenericModels/naivedecoder.jl" begin
        include("models/generic/naivedecoder.jl")
    end
end
