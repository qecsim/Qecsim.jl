using Test

using Qecsim.GenericModels
using Qecsim.Model
using Qecsim.BasicModels:FiveQubitCode
using Qecsim.PauliTools:to_pauli
using Random:MersenneTwister

# utility test function
function test_simpleerrormodel(
    error_model::SimpleErrorModel,
    p::Float64,
    expected_label::String,
    expected_probability_distribution::NTuple{4,Real},
    excluded_paulis::String=""
)
    # label
    @test label(error_model) == expected_label
    # probability_distribution
    @test probability_distribution(error_model, p) == expected_probability_distribution
    # generate
    code = FiveQubitCode()
    error = generate(error_model, code, p)
    error_pauli = to_pauli(error)
    @test length(error_pauli) == nkd(code)[1]
    for pauli in excluded_paulis
        @test count(==(pauli), error_pauli) == 0
    end
    # generate (seeded)
    error1 = generate(error_model, code, p, MersenneTwister(13))
    error2 = generate(error_model, code, p, MersenneTwister(13))
    @test error1 == error2
end

@testset "BitFlipErrorModel" begin
    error_model = BitFlipErrorModel()
    p = 0.1
    test_simpleerrormodel(error_model, p, "Bit-flip", (1 - p, p, 0, 0), "YZ")
end

@testset "BitPhaseFlipErrorModel" begin
    error_model = BitPhaseFlipErrorModel()
    p = 0.1
    test_simpleerrormodel(error_model, p, "Bit-phase-flip", (1 - p, 0, p, 0), "XZ")
end

@testset "DepolarizingErrorModel" begin
    error_model = DepolarizingErrorModel()
    p = 0.1
    test_simpleerrormodel(error_model, p, "Depolarizing", (1 - p, p / 3, p / 3, p / 3))
end

@testset "PhaseFlipErrorModel" begin
    error_model = PhaseFlipErrorModel()
    p = 0.1
    test_simpleerrormodel(error_model, p, "Phase-flip", (1 - p, 0, 0, p), "XY")
end
