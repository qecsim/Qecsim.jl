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
    # simple run
    code = FiveQubitCode()
    error_model = DepolarizingErrorModel()
    decoder = NaiveDecoder()
    p = 0.25
    max_runs = 1000
    data = qec_run(code, error_model, decoder, p; max_runs=max_runs)
    JSON.print(data, 4)
    expected_keys = Set([:code, :n_k_d, :time_steps, :error_model, :decoder,
        :error_probability, :measurement_error_probability, :n_run, :n_success, :n_fail,
        :n_logical_commutations, :custom_totals, :error_weight_total, :error_weight_pvar,
        :logical_failure_rate, :physical_error_rate, :wall_time])
    @test keys(data) == expected_keys
    @test data[:n_run] == max_runs
    @test data[:n_success] + data[:n_fail] == data[:n_run]
    @test data[:n_success] >= 0 && data[:n_fail] >= 0
    @test data[:n_fail] <= sum(data[:n_logical_commutations])
    @test data[:logical_failure_rate] == data[:n_fail] / data[:n_run]
    # physical_error_rate
    p_rate = data[:physical_error_rate]
    p_rate_std = sqrt(data[:error_weight_pvar] / (data[:n_k_d][1] ^ 2))
    @test p_rate - p_rate_std < p < p_rate + p_rate_std
    # run count
    data = qec_run(code, error_model, decoder, p)
    @test data[:n_run] == 1
    data = qec_run(code, error_model, decoder, p; max_runs=10)
    @test data[:n_run] == 10
    data = qec_run(code, error_model, decoder, p; max_failures=2)
    @test data[:n_fail] == 2
    data = qec_run(code, error_model, decoder, p; max_runs=10, max_failures=3)
    @test ((data[:n_run] == 10 && data[:n_fail] <= 3)
        || (data[:n_run] <= 10 && data[:n_fail] == 3))
    # seeded run
    data1 = qec_run(code, error_model, decoder, p, 13; max_runs=max_runs)
    data2 = qec_run(code, error_model, decoder, p, 13; max_runs=max_runs)
    delete!(data1, :wall_time)
    delete!(data2, :wall_time)
    @test data1 == data2

    #TODO: implement DecodeResult handling with parameterized custom_values numeric vector
    #      (see https://docs.julialang.org/en/v1/manual/performance-tips/#Type-declarations)
    #TODO: qec_run/run_once: support custom_values, custom_totals
    #TODO: qec_run docs
    #TODO: qec_merge
    #TODO: CLI/file versions of qec_run and qec_merge
    #TODO: test all major methods for type stability
end
