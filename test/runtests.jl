using SafeTestsets

@safetestset "Qecsim.jl" begin
    @safetestset "poc.jl" begin include("poc.jl") end
    @safetestset "PauliTools.jl" begin include("PauliTools.jl") end
    @safetestset "Model.jl" begin include("Model.jl") end
end
