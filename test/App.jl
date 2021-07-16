import Qecsim.Model

using Test
using Qecsim.App
using Qecsim.BasicModels: FiveQubitCode
using Qecsim.GenericModels: DepolarizingErrorModel, NaiveDecoder
using Qecsim.Model: ErrorModel, Decoder, DecodeResult
using Qecsim.PauliTools: to_bsf
using JSON
using Random: MersenneTwister

# test stub error model that generates a fixed given error
struct _FixedErrorModel <: ErrorModel
    error::BitVector
end
Model.label(::_FixedErrorModel) = "fixed"
Model.generate(error_model::_FixedErrorModel, x...) = error_model.error

# test stub decoder that decodes to a fixed given recovery
struct _FixedDecoder <: Decoder
    result::DecodeResult
end
Model.label(::_FixedDecoder) = "fixed"
Model.decode(decoder::_FixedDecoder, x...; kwargs...) = decoder.result

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
    # warn: RECOVERY DOES NOT RETURN TO CODESPACE
    error_model = _FixedErrorModel(to_bsf("IIIII"))
    decoder = _FixedDecoder(DecodeResult(to_bsf("IIIIX")))
    @test_logs (:warn, "RECOVERY DOES NOT RETURN TO CODESPACE") #=
        =# data = qec_run_once(code, error_model, decoder, p)
    @test !data[:success]
end

@testset "qec_run" begin
    code = FiveQubitCode()
    error_model = DepolarizingErrorModel()
    decoder = NaiveDecoder()
    p = 0.25
    data = qec_run(code, error_model, decoder, p; max_runs=1000)
    println("json=$(JSON.print(data, 4))")
end
