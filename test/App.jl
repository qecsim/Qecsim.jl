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

@testset "RunResult" begin
    r0 = RunResult(false, BitVector([0, 1]), 3)  # reference
    r1 = RunResult(false, BitVector([0, 1]), 3)  # same
    r2 = RunResult(true, BitVector([0, 1]), 3)   # different success
    r3 = RunResult(false, BitVector([0, 0]), 3)  # different logical_commutations
    r4 = RunResult(false, BitVector([0, 1]), 4)  # different error_weight
    @test r1 == r0 && hash(r1) == hash(r0)
    @test r2 != r0 && hash(r2) != hash(r0)
    @test r3 != r0 && hash(r3) != hash(r0)
    @test r4 != r0 && hash(r4) != hash(r0)
end

@testset "qec_run_once" begin
    # simple run
    code = FiveQubitCode()
    error_model = DepolarizingErrorModel()
    decoder = NaiveDecoder()
    p = 0.1
    data = qec_run_once(code, error_model, decoder, p)
    @test !data.success || !any(data.logical_commutations)  # success -> zero commutator
    # seeded run
    data1 = qec_run_once(code, error_model, decoder, p, MersenneTwister(13))
    data2 = qec_run_once(code, error_model, decoder, p, MersenneTwister(13))
    @test data1 == data2
    # warn: RECOVERY DOES NOT RETURN TO CODESPACE
    error_model = _FixedErrorModel(to_bsf("IIIII"))
    decoder = _FixedDecoder(DecodeResult(to_bsf("IIIIX")))
    @test_logs (:warn, "RECOVERY DOES NOT RETURN TO CODESPACE") #=
        =# data = qec_run_once(code, error_model, decoder, p)
    @test !data.success
end

@testset "qec_run" begin
    code = FiveQubitCode()
    error_model = DepolarizingErrorModel()
    decoder = NaiveDecoder()
    p = 0.25
    data = qec_run(code, error_model, decoder, p; max_runs=1000)
    println("json=$(JSON.print(data, 4))")
    #TODO: qec_run tests
    #TODO: implement DecodeResult handling with parameterized custom_values numeric vector
    #      (see https://docs.julialang.org/en/v1/manual/performance-tips/#Type-declarations)
    #TODO: qec_run/run_once: support custom_values, custom_totals
    #TODO: qec_run docs
    #TODO: qec_merge
    #TODO: CLI/file versions of qec_run and qec_merge
    #TODO: test all major methods for type stability
end
