using Test
using Qecsim.App
using Qecsim.BasicModels: FiveQubitCode
using Qecsim.GenericModels: DepolarizingErrorModel, NaiveDecoder
using Random: MersenneTwister

@testset "qec_run_once" begin
    # simple run
    code = FiveQubitCode()
    error_model = DepolarizingErrorModel()
    decoder = NaiveDecoder()
    p = 0.1
    data = qec_run_once(code, error_model, decoder, p)
    data_kt = Dict(
        :error_weight => Int,
        :logical_commutations => AbstractVector{Bool},
        :success => Bool,
    )
    @test keys(data) == keys(data_kt)  # test data keys
    @test all([typeof(v) <: data_kt[k] for (k, v) in data])  # test data values types
    @test !data[:success] || !any(data[:logical_commutations])  # success -> zero commutator
    # seeded run
    data1 = qec_run_once(code, error_model, decoder, p, MersenneTwister(13))
    data2 = qec_run_once(code, error_model, decoder, p, MersenneTwister(13))
    @test data1 == data2
    #TODO: test for @warn
end
