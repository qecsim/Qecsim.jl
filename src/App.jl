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
export RunResult, qec_run_once, qec_run

"""
    qec_run_once(code::StabilizerCode, error_model::ErrorModel, decoder::Decoder,
                 p::Float64, rng::AbstractRNG=GLOBAL_RNG)

Execute a stabilizer code error-decode-recovery (ideal) simulation and return run result.

TODO: complete doc
"""
function qec_run_once(
    code::StabilizerCode,
    error_model::ErrorModel,
    decoder::Decoder,
    p::Float64,
    rng::AbstractRNG=GLOBAL_RNG
)
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
    # apply overrides
    resolved_success = isnothing(success) ? s_commutes && l_commutes : success
    resolved_l_comms = isnothing(logical_commutations) ? l_comms : logical_commutations
    return (resolved_success, resolved_l_comms)
end

"""
    RunResult(success::Bool, error_weight::Int,
              logical_commutations::Union{Nothing,AbstractVector{Bool}}
              custom_values::Union{Nothing,Vector{<:Real}}
              )

Construct run result.

TODO: complete doc
"""
struct RunResult{T<:Union{Nothing,Vector{<:Real}}}
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


"""
    qec_run(code::StabilizerCode, error_model::ErrorModel, decoder::Decoder,
            p::Float64, random_seed;
            max_runs::Union{Int,Nothing}=nothing, max_failures::Union{Int,Nothing}=nothing)

Execute stabilizer code error-decode-recovery (ideal) simulations many times and return
aggregated run data.

TODO: complete doc
"""
function qec_run(
    code::StabilizerCode,
    error_model::ErrorModel,
    decoder::Decoder,
    p::Float64,
    random_seed=nothing;
    max_runs::Union{Int,Nothing}=nothing,
    max_failures::Union{Int,Nothing}=nothing,
)
    # derived defaults
    max_runs = isnothing(max_runs) && isnothing(max_failures) ? 1 : max_runs

    @info "qec_run: starting" code = code error_model = error_model decoder = decoder p = p #=
        =# random_seed = random_seed max_runs = max_runs max_failures = max_failures
    wall_time_start = time_ns()

    rng = MersenneTwister(random_seed)
    @info "qec_run: rng=$(rng)"

    # run counters
    n_run = n_success = 0
    error_weights = Int[]
    if !isnothing(max_runs) sizehint!(error_weights, max_runs) end
    n_logical_commutations::Union{Nothing,Vector{Int}} = nothing  # promote to Int[]
    custom_totals = nothing

    # do runs
    while ((isnothing(max_runs) || n_run < max_runs)
        && (isnothing(max_failures) || (n_run - n_success) < max_failures)
    )
        data = qec_run_once(code, error_model, decoder, p, rng)
        n_run += 1
        n_success += data.success ? 1 : 0
        push!(error_weights, data.error_weight)
        if n_run == 1  # initialize vector sums
            n_logical_commutations = _null_vec_copy(data.logical_commutations)
        else  # update vector sums
            _null_vec_add!(n_logical_commutations, data.logical_commutations)
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

end
