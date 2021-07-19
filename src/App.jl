"""
Functions to run quantum error correction simulations and merge output data.
"""
module App

# imports
using ..Model
using ..PauliTools: bsp, pack, weight
using Random: AbstractRNG, GLOBAL_RNG, MersenneTwister
using Statistics: var

# exports
export RunResult, qec_run_once, qec_run

"""
    qec_run_once(code::StabilizerCode, error_model::ErrorModel, decoder::Decoder,
        p::Float64, rng::AbstractRNG=GLOBAL_RNG)

Run a stabilizer code error-decode-recovery (ideal) simulation and return run data.
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
    recovered = xor.(error, result.recovery)
    @debug "qec_run_once: recovered=$(recovered)" recovered
    s_commutations = bsp(stabilizers(code), recovered)
    l_commutations = bsp(logicals(code), recovered)
    @debug "qec_run_once: stablizer_commutations=$(s_commutations)" s_commutations
    @debug "qec_run_once: logical_commutations=$(l_commutations)" l_commutations
    s_commutes = !any(s_commutations)
    l_commutes = !any(l_commutations)
    if !s_commutes
        log_data = Dict(  # pack to concise string format
            :error => pack(error),
            :recovery => pack(result.recovery),
        )
        @warn "RECOVERY DOES NOT RETURN TO CODESPACE" log_data...
    end
    success = s_commutes && l_commutes
    @debug "qec_run_once: success=$(success)"
    return RunResult(success, l_commutations, weight(error))
end
struct RunResult
    success::Bool
    logical_commutations::Vector{Bool}
    error_weight::Int
end
# equality methods TODO: write equality macro
function Base.hash(a::RunResult, h::UInt)
    h = hash(typeof(a), h)
    h = hash(a.success, h)
    h = hash(a.logical_commutations, h)
    return hash(a.error_weight, h)
end
function Base.:(==)(a::RunResult, b::RunResult)
    return (a.success == b.success
        && a.logical_commutations == b.logical_commutations
        && a.error_weight == b.error_weight)
end
function Base.isequal(a::RunResult, b::RunResult)
    return (isequal(a.success, b.success)
        && isequal(a.logical_commutations, b.logical_commutations)
        && isequal(a.error_weight, b.error_weight))
end


function qec_run(
    code::StabilizerCode,
    error_model::ErrorModel,
    decoder::Decoder,
    p::Float64,
    random_seed=nothing;
    max_runs::Union{Int, Nothing}=nothing,
    max_failures::Union{Int, Nothing}=nothing,
)
    # derived defaults
    max_runs = isnothing(max_runs) && isnothing(max_failures) ? 1 : max_runs

    @info "qec_run: starting" code=code error_model=error_model decoder=decoder p=p #=
        =# random_seed=random_seed max_runs=max_runs max_failures=max_failures
    wall_time_start = time_ns()

    rng = MersenneTwister(random_seed)
    @info "qec_run: rng=$(rng)"

    # run counters
    n_run = n_success = 0
    error_weights = Int[]
    if !isnothing(max_runs) sizehint!(error_weights, max_runs) end
    n_logical_commutations::Union{Nothing, Vector{Int}} = nothing  # promote to Int[]
    custom_totals = nothing

    # do runs
    while ((isnothing(max_runs) || n_run < max_runs)
           && (isnothing(max_failures) || (n_run - n_success) < max_failures))
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
