using Test
using Qecsim.App
using Qecsim.BasicModels:FiveQubitCode
using Qecsim.GenericModels:BitFlipErrorModel, DepolarizingErrorModel, NaiveDecoder
using Qecsim.Model:Model, ErrorModel, Decoder, DecodeResult, logical_xs, logical_zs
using Qecsim.PauliTools:to_bsf
using JSON
using Random:MersenneTwister

# TODO: qec_read/qec_write doc
# TODO: profiling / type-stability checks


# test stub error model that generates a fixed given error
struct _FixedErrorModel <: ErrorModel
    error::BitVector
end
Model.label(::_FixedErrorModel) = "fixed"
Model.generate(error_model::_FixedErrorModel, x...) = error_model.error
# test stub decoder that decodes to a fixed given decode result
struct _FixedDecoder <: Decoder
    result::DecodeResult
end
Model.label(::_FixedDecoder) = "fixed"
Model.decode(decoder::_FixedDecoder, x...; kwargs...) = decoder.result
# test stub decoder that decodes to a cycle of given decode results
struct _CycleDecoder <: Decoder
    result_iter  # Iterator{DecodeResult}
    _CycleDecoder(results) = new(Iterators.Stateful(Iterators.Cycle(results)))
end
Model.label(::_CycleDecoder) = "cycle"
Model.decode(decoder::_CycleDecoder, x...; kwargs...) = popfirst!(decoder.result_iter)


# test stub code using duck-typing
struct _DuckCode end
Model.label(::_DuckCode) = "duck"
Model.stabilizers(::_DuckCode) = to_bsf(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"])
Model.logical_xs(::_DuckCode) = to_bsf(["XXXXX"])
Model.logical_zs(::_DuckCode) = to_bsf(["ZZZZZ"])
Model.nkd(::_DuckCode) = (5, 1, 3)
# test stub error model using duck-typing
struct _DuckErrorModel end
Model.label(::_DuckErrorModel) = "duck"
Model.generate(::_DuckErrorModel, x...) = to_bsf("IIIII")
# test stub decoder using duck-typing
struct _DuckDecoder end
Model.label(::_DuckDecoder) = "duck"
Model.decode(::_DuckDecoder, x...; kwargs...) = DecodeResult(recovery=to_bsf("IIIII"))


@testset "RunResult" begin
    r0 = RunResult(false, 3, BitVector([0, 1]), nothing)  # reference
    r1 = RunResult(false, 3, BitVector([0, 1]), nothing)  # same
    r2 = RunResult(true, 3, BitVector([0, 1]), nothing)   # different success
    r3 = RunResult(false, 3, BitVector([0, 0]), nothing)  # different logical_commutations
    r4 = RunResult(false, 4, BitVector([0, 1]), nothing)  # different error_weight
    @test r1 == r0 && hash(r1) == hash(r0)
    @test r2 != r0 && hash(r2) != hash(r0)
    @test r3 != r0 && hash(r3) != hash(r0)
    @test r4 != r0 && hash(r4) != hash(r0)
    @test isequal(r1, r0) && hash(r1) == hash(r0)
    @test !isequal(r2, r0) && hash(r2) != hash(r0)
    @test !isequal(r3, r0) && hash(r3) != hash(r0)
    @test !isequal(r4, r0) && hash(r4) != hash(r0)
end

@testset "qec_run_once" begin
    # simulation parameters
    code = FiveQubitCode()
    error_model = DepolarizingErrorModel()
    decoder = NaiveDecoder()
    # float p
    data = qec_run_once(code, error_model, decoder, 0.1)
    @test !data.success || !any(data.logical_commutations)  # success -> zero commutator
    # rational p
    data = qec_run_once(code, error_model, decoder, 1 // 10)
    @test !data.success || !any(data.logical_commutations)  # success -> zero commutator
    # seeded run
    data1 = qec_run_once(code, error_model, decoder, 0.1, MersenneTwister(13))
    data2 = qec_run_once(code, error_model, decoder, 0.1, MersenneTwister(13))
    @test data1 == data2
    # warn: RECOVERY DOES NOT RETURN TO CODESPACE
    error_model = _FixedErrorModel(to_bsf("IIIII"))
    decoder = _FixedDecoder(DecodeResult(recovery=to_bsf("IIIIX")))
    @test_logs (:warn, "RECOVERY DOES NOT RETURN TO CODESPACE") #=
        =# data = qec_run_once(code, error_model, decoder, 0.1)
    @test !data.success
end

@testset "qec_run_once-duck-typing" begin
    # duck-type simulation
    data = qec_run_once(_DuckCode(), _DuckErrorModel(), _DuckDecoder(), 0.1)
    @test data.success
    @test !any(data.logical_commutations)
    @test data.error_weight == 0
    @test isnothing(data.custom_values)
end

@testset "qec_run_once-override" begin
    # common simulation parameters
    code = FiveQubitCode()
    identity = to_bsf("IIIII")
    error_model = _FixedErrorModel(identity)
    p = 0.1
    # tests with and without overrides
    for (decode_result, expected) in [
        # identity recovery
        (DecodeResult(recovery=identity),
            RunResult(true, 0, BitVector([0, 0]), nothing)),
        # logical_x recovery
        (DecodeResult(recovery=logical_xs(code)[1,:]),
            RunResult(false, 0, BitVector([0, 1]), nothing)),
        # logical_z recovery
        (DecodeResult(recovery=logical_zs(code)[1,:]),
            RunResult(false, 0, BitVector([1, 0]), nothing)),
        # identity but override success=false
        (DecodeResult(success=false, recovery=identity),
            RunResult(false, 0, BitVector([0, 0]), nothing)),
        # identity but override logical_commutations=[1, 1]
        (DecodeResult(recovery=identity, logical_commutations=[1, 1]),
            RunResult(true, 0, BitVector([1, 1]), nothing)),
        # identity but override success=false and logical_commutations=[1, 1]
        (DecodeResult(success=false, recovery=identity, logical_commutations=[1, 1]),
            RunResult(false, 0, BitVector([1, 1]), nothing)),
        # no-recovery but override success=false
        (DecodeResult(success=false),
            RunResult(false, 0, nothing, nothing)),
        # no-recovery but override success=false and logical_commutations=[1, 1]
        (DecodeResult(success=false, logical_commutations=[1, 1]),
            RunResult(false, 0, BitVector([1, 1]), nothing)),
        # no-recovery but override success=false and custom_values=[1, 2]
        (DecodeResult(success=false, custom_values=[1, 2]),
            RunResult(false, 0, nothing, [1, 2])),
        # no-recovery but override success=false and custom_values=[1., 2.]
        (DecodeResult(success=false, custom_values=[1., 2.]),
            RunResult(false, 0, nothing, [1., 2.])),
        # no-recovery but override success=false and custom_values=[1//2, 1//3]
        (DecodeResult(success=false, custom_values=[1 // 2, 1 // 3]),
            RunResult(false, 0, nothing, [1 // 2, 1 // 3])),
    ]
        decoder = _FixedDecoder(decode_result)
        data = qec_run_once(code, error_model, decoder, p)
        @test data == expected
    end
end

@testset "qec_run" begin
    # simulation parameters
    code = FiveQubitCode()
    error_model = DepolarizingErrorModel()
    decoder = NaiveDecoder()
    max_runs = 1000
    for p in [0.25, 1 // 4]  # test with float and rational p
        data = qec_run(code, error_model, decoder, p; max_runs=max_runs)
        # JSON.print(data, 4)
        expected_keys = Set([:code, :n_k_d, :time_steps, :error_model, :decoder,
            :error_probability, :measurement_error_probability, :n_run, :n_success, :n_fail,
            :n_logical_commutations, :custom_totals, :error_weight_total,
            :error_weight_pvar, :logical_failure_rate, :physical_error_rate, :wall_time])
        @test keys(data) == expected_keys
        @test data[:error_probability] == p
        @test typeof(data[:error_probability]) == typeof(p)
        @test data[:n_run] == max_runs
        @test data[:n_success] + data[:n_fail] == data[:n_run]
        @test data[:n_success] >= 0 && data[:n_fail] >= 0
        @test data[:n_fail] <= sum(data[:n_logical_commutations])
        @test data[:logical_failure_rate] == data[:n_fail] / data[:n_run]
        # physical_error_rate
        p_rate = data[:physical_error_rate]
        p_rate_std = sqrt(data[:error_weight_pvar] / (data[:n_k_d][1]^2))
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
    end
    # seeded run
    data1 = qec_run(code, error_model, decoder, 0.1, 13; max_runs=max_runs)
    data2 = qec_run(code, error_model, decoder, 0.1, 13; max_runs=max_runs)
    delete!(data1, :wall_time)
    delete!(data2, :wall_time)
    @test data1 == data2
end

@testset "qec_run-duck-typing" begin
    # duck-type simulation
    data = qec_run(_DuckCode(), _DuckErrorModel(), _DuckDecoder(), 0.1; max_runs=10)
    @test data[:n_success] == 10
    @test all(data[:n_logical_commutations] .== 0)
    @test data[:error_weight_total] == 0
    @test isnothing(data[:custom_totals])
end

@testset "qec_run-override" begin
    # common simulation parameters
    code = FiveQubitCode()
    identity = to_bsf("IIIII")
    error_model = _FixedErrorModel(identity)
    p = 0.1
    # tests with and without overrides
    for (decoder, max_runs, expected) in [
        # identity recovery
        (_FixedDecoder(DecodeResult(recovery=identity)), 1,
            Dict(:n_success => 1, :n_fail => 0, :n_logical_commutations => [0, 0])),
        # logical_x recovery
        (_FixedDecoder(DecodeResult(recovery=logical_xs(code)[1,:])), 2,
            Dict(:n_success => 0, :n_fail => 2, :n_logical_commutations => [0, 2])),
        # logical_z recovery
        (_FixedDecoder(DecodeResult(recovery=logical_zs(code)[1,:])), 3,
            Dict(:n_success => 0, :n_fail => 3, :n_logical_commutations => [3, 0])),
        # identity but override success=false
        (_FixedDecoder(DecodeResult(success=false, recovery=identity)), 4,
            Dict(:n_success => 0, :n_fail => 4, :n_logical_commutations => [0, 0])),
        # identity but override logical_commutations=[1, 1]
        (_FixedDecoder(DecodeResult(recovery=identity, logical_commutations=[1, 1])), 5,
            Dict(:n_success => 5, :n_fail => 0, :n_logical_commutations => [5, 5])),
        # identity but override success=false and logical_commutations=[1, 1]
        (_FixedDecoder(DecodeResult(success=false, recovery=identity,
            logical_commutations=[1, 1])), 6,
            Dict(:n_success => 0, :n_fail => 6, :n_logical_commutations => [6, 6])),
        # no-recovery but override success=false
        (_FixedDecoder(DecodeResult(success=false)), 7,
            Dict(:n_success => 0, :n_fail => 7, :n_logical_commutations => nothing)),
        # no-recovery but override success=false and logical_commutations=[1, 1]
        (_FixedDecoder(DecodeResult(success=false, logical_commutations=[1, 1])), 8,
            Dict(:n_success => 0, :n_fail => 8, :n_logical_commutations => [8, 8])),
        # no-recovery but override success=false and custom_values=[1, 1]
        (_FixedDecoder(DecodeResult(success=false, custom_values=[1, 1])), 9,
            Dict(:n_success => 0, :n_fail => 9, :custom_totals => [9, 9])),
        # no-recovery but override success=false and custom_values=[1., 1.]
        (_FixedDecoder(DecodeResult(success=false, custom_values=[1., 1.])), 10,
            Dict(:n_success => 0, :n_fail => 10, :custom_totals => [10., 10.])),
        # no-recovery but override success=false and custom_values=[1//2, 1//3]
        (_FixedDecoder(DecodeResult(success=false, custom_values=[1 // 2, 1 // 3])), 11,
            Dict(:n_success => 0, :n_fail => 11, :custom_totals => [11 // 2, 11 // 3])),
        # no-recovery but override success=false/true and custom_values=[1, 1]/[2, 3]
        (_CycleDecoder([DecodeResult(success=false, custom_values=[1, 1]),
            DecodeResult(success=true, custom_values=[2, 3])]), 12,
            Dict(:n_success => 6, :n_fail => 6, :custom_totals => [18, 24])),
    ]
        data = qec_run(code, error_model, decoder, p; max_runs=max_runs)
        @test issubset(expected, data)
    end
    # test that inconsistent overrides fail
    for (decode_results, expected) in [
        ([DecodeResult(success=true, logical_commutations=nothing),
            DecodeResult(success=true, logical_commutations=[1, 0])], MethodError),
        ([DecodeResult(success=true, logical_commutations=[1, 1, 0]),
            DecodeResult(success=true, logical_commutations=[1, 0])], DimensionMismatch),
        ([DecodeResult(success=true, custom_values=[1, 0]),
            DecodeResult(success=true, custom_values=[1.1, 0.1])], InexactError),
    ]
        decoder = _CycleDecoder(decode_results)
        @test_throws expected qec_run(code, error_model, decoder, p; max_runs=5)
    end
end

@testset "qec_merge-similar" begin
    # similar simulation parameters
    code = FiveQubitCode()
    error_model = DepolarizingErrorModel()
    decoder = NaiveDecoder()
    p = 0.20
    max_runs = 10
    data1 = qec_run(code, error_model, decoder, p; max_runs=max_runs)
    # JSON.print(data1, 4)
    data2 = qec_run(code, error_model, decoder, p; max_runs=max_runs)
    # JSON.print(data2, 4)
    # data1 and data2 are dicts of same length
    @test length(data1) == length(data2)
    merged_data_list = qec_merge(data1, data2)
    # JSON.print(merged_data_list, 4)
    # merged_data_list is a vector with one entry
    @test length(merged_data_list) == 1
    merged_data = merged_data_list[1]
    # merged_data drops the error_weight_pvar key
    @test length(merged_data) == length(data1) - 1
    # merged_data preserves p and sums max_runs
    @test merged_data[:error_probability] == p
    @test merged_data[:n_run] == max_runs * 2
    # all value types of merged_data match corresponding value types of data1
    @test all(typeof(merged_data[k]) == typeof(data1[k]) for k in keys(merged_data))
end

@testset "qec_merge-distinct" begin
    # distinct simulation parameters
    code = FiveQubitCode()
    error_model = DepolarizingErrorModel()
    decoder = NaiveDecoder()
    p1, p2 = 0.1, 0.2
    max_runs = 10
    data1 = qec_run(code, error_model, decoder, p1; max_runs=max_runs)
    data2 = qec_run(code, error_model, decoder, p2; max_runs=max_runs)
    # data1 and data2 are dicts of same length
    @test length(data1) == length(data2)
    merged_data_list = qec_merge(data1, data2)
    # merged_data_list is a vector with two entries
    @test length(merged_data_list) == 2
    merged_data1 = merged_data_list[1]
    merged_data2 = merged_data_list[2]
    # merged_data drops the error_weight_pvar key
    @test length(merged_data1) == length(data1) - 1
    @test length(merged_data2) == length(data2) - 1
    # merged_data preserves all values except :error_weight_pvar
    delete!(data1, :error_weight_pvar)
    @test merged_data1 == data1
    delete!(data2, :error_weight_pvar)
    @test merged_data2 == data2
end

@testset "qec_merge-zero" begin
    # no data to merge
    merged_data_list = qec_merge()
    # merged_data_list is a vector with no entries
    @test length(merged_data_list) == 0
end

@testset "qec_merge-one" begin
    # single simulation
    data = qec_run(FiveQubitCode(), BitFlipErrorModel(), NaiveDecoder(), 0.1; max_runs=10)
    merged_data_list = qec_merge(data)
    # merged_data_list is a vector with one entry
    @test length(merged_data_list) == 1
    merged_data = merged_data_list[1]
    # merged_data drops the error_weight_pvar key
    @test length(merged_data) == length(data) - 1
    # merged_data preserves all values except :error_weight_pvar
    delete!(data, :error_weight_pvar)
    @test merged_data == data
end

@testset "qec_write_read" begin
    # multiple simulation parameters
    code = FiveQubitCode()
    error_model = DepolarizingErrorModel()
    decoder = NaiveDecoder()
    p1, p2 = 0.08, 0.10
    max_runs = 10
    data = []
    push!(data, qec_run(code, error_model, decoder, p1; max_runs=max_runs))
    push!(data, qec_run(code, error_model, decoder, p2; max_runs=max_runs))
    # JSON.print(data, 4)
    # write data
    filename = tempname()
    qec_write(filename, data...)
    # read data
    data_read = qec_read(filename)
    # JSON.print(data_read, 4)
    @test data_read == data
end

@testset "qec_write-existing-file" begin
    # single simulation
    data = qec_run(FiveQubitCode(), BitFlipErrorModel(), NaiveDecoder(), 0.1; max_runs=10)
    # write data
    filename = tempname()
    touch(filename)
    # refuse to overwrite existing file
    @test_throws ErrorException qec_write(filename, data...)
end

@testset "qec_write_read-empty-data" begin
    # write empty data
    filename = tempname()
    qec_write(filename)
    # read empty data
    data_read = qec_read(filename)
    @test data_read == Dict{Symbol,Any}[]
end

@testset "qec_read-invalid-data" begin
    # write invalid data
    filename = tempname()
    touch(filename)
    # fail read invalid data
    @test_throws Exception qec_read(filename)
    # write invalid data
    filename = tempname()
    open(io -> JSON.print(io, Dict(1 => '1')), filename, "w")
    # fail read invalid data
    @test_throws Exception qec_read(filename)
end
