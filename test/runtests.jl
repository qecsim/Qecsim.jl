using SafeTestsets

@safetestset "poc" begin include("poc.jl") end
@safetestset "PauliTools" begin include("PauliTools/PauliTools.jl") end
