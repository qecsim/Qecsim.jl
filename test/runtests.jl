using SafeTestsets

@safetestset "Qecsim.jl" begin
    @safetestset "poc.jl" begin include("poc.jl") end
    @safetestset "PauliTools.jl" begin include("PauliTools/PauliTools.jl") end
end
