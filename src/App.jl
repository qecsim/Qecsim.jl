"""
Functions to run quantum error correction simulations and merge output data.
"""
module App

# imports
using ..Model
using ..PauliTools:bsp, pack, weight
using Random:AbstractRNG, GLOBAL_RNG, MersenneTwister
using Statistics:var

# exports
export RunResult, qec_merge, qec_run_once, qec_run

@doc raw"""
    qec_run_once(code, error_model, decoder, p::Real, rng::AbstractRNG=GLOBAL_RNG)
        -> RunResult

Execute a stabilizer code error-decode-recovery (ideal) simulation and return run result.

The parameters `code`, `error_model` and `decoder` should be concrete subtypes or duck-typed
implementations of [`StabilizerCode`](@ref), [`ErrorModel`](@ref) and [`Decoder`](@ref),
respectively.

The simulation algorithm is as follows:

1. ``S ←`` `stabilizers(code)`
2. ``L ←`` `logicals(code)`
3. ``e ←`` `generate(error_model, code, p, rng)`
4. ``y ← S ⊙ e``
5. `decode_result` ``←`` `decode(decoder, code,` ``y```; kwargs...)`
6. ``r ←`` `decode_result.recovery`
7. sanity check: ``S ⊙ (r ⊕ e) = 0``
8. `logical_commutations` ``← L ⊙ (r ⊕ e)``
8. `success` ``← L ⊙ (r ⊕ e) = 0``

where ``⊕`` denotes element-wise exclusive-or, and ``⊙`` is defined in
[`PauliTools.bsp`](@ref bsp).

The `kwargs` passed to [`decode`](@ref) include `error_model`, `p` and `error`; most
decoders will ignore these parameters. The [`decode`](@ref) method returns a
[`DecodeResult`](@ref). If `decode_result.success` and/or
`decode_result.logical_commutations` are specified, they override the values of `success`
and `logical_commutations`, irrespective of whether `decode_result.recovery` is specified or
not. The value `decode_result.custom_values` is passed through in the run result.

See also [`RunResult`](@ref).

# Examples
```jldoctest
julia> using Qecsim.BasicModels, Qecsim.GenericModels, Random

julia> rng = MersenneTwister(6);  # use random seed for reproducible result

julia> qec_run_once(FiveQubitCode(), DepolarizingErrorModel(), NaiveDecoder(), 0.2, rng)
RunResult{Nothing}(false, 2, Bool[1, 0], nothing)
```
"""
function qec_run_once(code, error_model, decoder, p::Real, rng::AbstractRNG=GLOBAL_RNG)
    error = generate(error_model, code, p, rng)
    @debug "qec_run_once: error=$(error)" error
    syndrome = bsp(stabilizers(code), error)
    @debug "qec_run_once: syndrome=$(syndrome)" syndrome
    ctx = Dict(:error_model => error_model, :p => p, :error => error)
    result = decode(decoder, code, syndrome; ctx...)
    @debug "qec_run_once: result=$(result)" result
    success, l_commutations = _resolve_decoding(result.recovery, result.success,
        result.logical_commutations, code, error)
    @debug "qec_run_once: success=$(success)"
    @debug "qec_run_once: logical_commutations=$(l_commutations)"
    return RunResult(success, weight(error), l_commutations, result.custom_values)
end

# resolve tuple (success, logical_commutations)
function _resolve_decoding(::Nothing, success, logical_commutations, code, error)
    # recovery is null, so just return overrides
    return (success, logical_commutations)
end
function _resolve_decoding(recovery, success, logical_commutations, code, error)
    # recovery is not null, so evaluate success and logical_commutations
    recovered = xor.(error, recovery)
    s_comms = bsp(stabilizers(code), recovered)
    l_comms = bsp(logicals(code), recovered)
    s_commutes = !any(s_comms)
    l_commutes = !any(l_comms)
    if !s_commutes
        @warn "RECOVERY DOES NOT RETURN TO CODESPACE" pack(error) pack(recovery)
    end
    # derived overrides
    dsuccess = isnothing(success) ? s_commutes && l_commutes : success
    dl_comms = isnothing(logical_commutations) ? l_comms : logical_commutations
    return (dsuccess, dl_comms)
end

"""
    RunResult(success::Bool, error_weight::Int,
              logical_commutations::Union{Nothing,BitVector}
              custom_values::Union{Nothing,Vector{<:Real}})

Construct a run result as returned by [`qec_run_once`](@ref).

# Examples
```jldoctest
julia> r = RunResult(false, 2, BitVector([0, 1]), [1.2, 3.1])
RunResult{Vector{Float64}}(false, 2, Bool[0, 1], [1.2, 3.1])

julia> r.success, r.error_weight, r.logical_commutations, r.custom_values
(false, 2, Bool[0, 1], [1.2, 3.1])
```
"""
struct RunResult{T <: Union{Nothing,Vector{<:Real}}}
    success::Bool
    error_weight::Int
    logical_commutations::Union{Nothing,BitVector}
    custom_values::T
end
# equality methods TODO: write equality macro
function Base.hash(a::RunResult, h::UInt)
    h = hash(typeof(a), h)
    h = hash(a.success, h)
    h = hash(a.error_weight, h)
    h = hash(a.logical_commutations, h)
    h = hash(a.custom_values, h)
    return h
end
function Base.:(==)(a::RunResult, b::RunResult)
    return (a.success == b.success
        && a.error_weight == b.error_weight
        && a.logical_commutations == b.logical_commutations
        && a.custom_values == b.custom_values)
end
function Base.isequal(a::RunResult, b::RunResult)
    return (isequal(a.success, b.success)
        && isequal(a.error_weight, b.error_weight)
        && isequal(a.logical_commutations, b.logical_commutations)
        && isequal(a.custom_values, b.custom_values))
end


@doc raw"""
    qec_run(code, error_model, decoder, p::Real, random_seed=nothing;
            max_runs::Union{Integer,Nothing}=nothing,
            max_failures::Union{Integer,Nothing}=nothing)
        -> Dict

Execute stabilizer code error-decode-recovery (ideal) simulations many times and return
aggregated run data, see [`qec_run_once`](@ref) for details of a single run.

The parameters `code`, `error_model` and `decoder` should be concrete subtypes or duck-typed
implementations of [`StabilizerCode`](@ref), [`ErrorModel`](@ref) and [`Decoder`](@ref),
respectively.

Simulations are run one or more times as determined by `max_runs` and `max_failures`. If
`max_runs` and/or `max_failures` are specified, stop after `max_runs` runs or `max_failures`
failures, whichever happens first. If neither is specified, stop after one run.

The returned aggregated data has the following format:

    Dict(
        :code => "5-qubit"              # label(code)
        :n_k_d => (5, 1, 3)             # nkd(code)
        :time_steps => 1                # always 1 for ideal simulations
        :error_model => "Depolarizing"  # label(error_model)
        :decoder => "Naive"             # label(decoder)
        :error_probability => 0.1       # p
        :measurement_error_probability => 0.0   # always 0.0 for ideal simulations
        :n_run => 100                   # count of runs
        :n_success => 92                # count of successful recoveries
        :n_fail => 8                    # count of failed recoveries
        :n_logical_commutations => [5, 6]   # count of logical_commutations
        :custom_totals => nothing       # sum of custom_values
        :error_weight_total => 55       # sum of error_weight over n_run runs
        :error_weight_pvar => 0.4075    # pvariance of error_weight over n_run runs
        :logical_failure_rate => 0.08   # n_fail / n_run
        :physical_error_rate => 0.11    # error_weight_total / n_k_d[1] / time_steps / n_run
        :wall_time => 0.00253906        # wall-time for run in fractional seconds
    )

# Examples
```jldoctest; filter = r":wall_time +=> \d*\.?\d*"
julia> using Qecsim.BasicModels, Qecsim.GenericModels

julia> seed = 7;

julia> data = qec_run(FiveQubitCode(), DepolarizingErrorModel(), NaiveDecoder(), 0.1, seed;
           max_runs=100);
┌ Info: qec_run: starting
│   code = Qecsim.BasicModels.BasicCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"], ["XXXXX"], ["ZZZZZ"], (5, 1, 3), "5-qubit")
│   error_model = Qecsim.GenericModels.DepolarizingErrorModel()
│   decoder = Qecsim.GenericModels.NaiveDecoder(10)
│   p = 0.1
│   random_seed = 7
│   max_runs = 100
└   max_failures = nothing
[ Info: qec_run: rng=MersenneTwister(7)
[ Info: qec_run: complete: data=Dict{Symbol, Any}(:error_weight_pvar => 0.4075000000000001, :time_steps => 1, :n_logical_commutations => [5, 6], :error_weight_total => 55, :wall_time => 0.002539058, :n_k_d => (5, 1, 3), :error_model => "Depolarizing", :physical_error_rate => 0.11, :measurement_error_probability => 0.0, :error_probability => 0.1, :n_success => 92, :logical_failure_rate => 0.08, :custom_totals => nothing, :code => "5-qubit", :decoder => "Naive", :n_fail => 8, :n_run => 100)

julia> data
Dict{Symbol, Any} with 17 entries:
  :error_weight_pvar             => 0.4075
  :time_steps                    => 1
  :n_logical_commutations        => [5, 6]
  :error_weight_total            => 55
  :wall_time                     => 0.00253906
  :n_k_d                         => (5, 1, 3)
  :error_model                   => "Depolarizing"
  :physical_error_rate           => 0.11
  :measurement_error_probability => 0.0
  :error_probability             => 0.1
  :n_success                     => 92
  :logical_failure_rate          => 0.08
  :custom_totals                 => nothing
  :code                          => "5-qubit"
  :decoder                       => "Naive"
  :n_fail                        => 8
  :n_run                         => 100
```
"""
function qec_run(code, error_model, decoder, p::Real, random_seed=nothing;
    max_runs::Union{Integer,Nothing}=nothing, max_failures::Union{Integer,Nothing}=nothing,
)
    # derived defaults
    dmax_runs = isnothing(max_runs) && isnothing(max_failures) ? 1 : max_runs

    @info "qec_run: starting" code = code error_model = error_model decoder = decoder #=
        =# p = p random_seed = random_seed max_runs = dmax_runs max_failures = max_failures
    wall_time_start = time_ns()

    rng = MersenneTwister(random_seed)
    @info "qec_run: rng=$(rng)"

    # run counters
    n_run = n_success = 0
    error_weights = Int[]
    if !isnothing(dmax_runs) sizehint!(error_weights, dmax_runs) end
    n_logical_commutations::Union{Nothing,Vector{Int}} = nothing  # promote to Int[]
    custom_totals = nothing

    # do runs
    while ((isnothing(dmax_runs) || n_run < dmax_runs)
        && (isnothing(max_failures) || (n_run - n_success) < max_failures)
    )
        data = qec_run_once(code, error_model, decoder, p, rng)
        n_run += 1
        n_success += data.success ? 1 : 0
        push!(error_weights, data.error_weight)
        if n_run == 1  # initialize vector sums
            n_logical_commutations = _null_vec_copy(data.logical_commutations)
            custom_totals = _null_vec_copy(data.custom_values)
        else  # update vector sums
            _null_vec_add!(n_logical_commutations, data.logical_commutations)
            _null_vec_add!(custom_totals, data.custom_values)
        end
    end

    # prepare data
    runs_data = Dict(
        :code => label(code),
        :n_k_d => nkd(code),
        :time_steps => 1,
        :error_model => label(error_model),
        :decoder => label(decoder),
        :error_probability => p,
        :measurement_error_probability => 0.0,
        :n_run => n_run,
        :n_success => n_success,
        :n_fail => n_run - n_success,
        :n_logical_commutations => n_logical_commutations,
        :custom_totals => custom_totals,
        :error_weight_total => sum(error_weights),
        :error_weight_pvar => var(error_weights; corrected=false),
        :logical_failure_rate => 0.0,  # set below
        :physical_error_rate => 0.0,  # set below
        :wall_time => 0.0,  # set below
    )
    # set rate statistics
    _rate_statistics!(runs_data)

    runs_data[:wall_time] = (time_ns() - wall_time_start) / 10^9
    @info "qec_run: complete: data=$(runs_data)"

    return runs_data
end

# return copy of val (or nothing if null)
_null_vec_copy(::Nothing) = nothing
_null_vec_copy(val::AbstractVector) = copy(val)
# update total adding val (or nothing if both null; MethodError if one and only one is null)
_null_vec_add!(::Nothing, ::Nothing) = nothing
_null_vec_add!(total::AbstractVector, val::AbstractVector) = total .+= val

# set :logical_failure_rate and :physical_error_rate as defined in qec_run
function _rate_statistics!(runs_data)
    # extract data
    time_steps = runs_data[:time_steps]
    n_run = runs_data[:n_run]
    n_fail = runs_data[:n_fail]
    error_weight_total = runs_data[:error_weight_total]
    nkd = runs_data[:n_k_d][1]
    # add rate statistics
    runs_data[:logical_failure_rate] = n_fail / n_run
    runs_data[:physical_error_rate] = error_weight_total / nkd[1] / time_steps / n_run
    return nothing
end

function qec_merge(data...)
    # define group keys, value keys and zero values
    grp_keys = (:code, :n_k_d, :error_model, :decoder, :error_probability, :time_steps,
                :measurement_error_probability)
    scalar_val_keys = (:n_run, :n_fail, :n_success, :error_weight_total, :wall_time)
    scalar_zero_vals = (0, 0, 0, 0, 0.0)
    vector_val_keys = (:n_logical_commutations, :custom_totals,)
    # map of groups to sums
    grps_to_scalar_sums = Dict()
    grps_to_vector_sums = Dict()
    # iterate through single list from given data lists
    for runs_data in data
        # define defaults, create new data with defaults overwritten by data
        # support for 0.10 and 0.15 files:
        defaults_0_16 = Dict(:time_steps => 1, :measurement_error_probability => 0.0)
        # support for pre-1.0b6 files:
        defaults_1_0b6 = Dict(:n_logical_commutations => nothing, :custom_totals => nothing)
        # runs_data = Dict(defaults_0_16..., defaults_1_0b6..., runs_data...)
        runs_data = merge(defaults_0_16, defaults_1_0b6, runs_data)
        # extract group from data
        grp_id = Tuple(runs_data[k] for k in grp_keys)
        # scalars: e.g. (10, 6, 4, 256, 10.34) extracted from data
        scalar_vals = Tuple(runs_data[k] for k in scalar_val_keys)
        scalar_sums = get(grps_to_scalar_sums, grp_id, scalar_zero_vals)  # get sums
        scalar_sums = Tuple(sum(x) for x in zip(scalar_sums, scalar_vals))  # update sums
        grps_to_scalar_sums[grp_id] = scalar_sums  # put sums
        # vectors: e.g. ([2, 5], [3, 8, 2], nothing) extracted from data
        vector_vals = Tuple(_null_vec_copy(runs_data[k]) for k in vector_val_keys)
        if haskey(grps_to_vector_sums, grp_id)  # sums already in map
            vector_sums = grps_to_vector_sums[grp_id]  # get sums
            vector_sums = Tuple(_null_vec_add!(s, v)  # update sums
                                for (s, v) in zip(vector_sums, vector_vals))
        else  # sums not in map yet
            vector_sums = vector_vals  # update sums
        end
        grps_to_vector_sums[grp_id] = vector_sums  # put sums
    end
    # flatten grps_to_scalar_sums and grps_to_vector_sums
    merged_data_list = [Dict(
        zip((grp_keys..., scalar_val_keys..., vector_val_keys...),
            (grp_id..., scalar_sums..., grps_to_vector_sums[grp_id]...))
        ) for (grp_id, scalar_sums) in grps_to_scalar_sums]
    # update rate statistics
    for runs_data in merged_data_list
        _rate_statistics!(runs_data)
    end
    return merged_data_list
end

end
