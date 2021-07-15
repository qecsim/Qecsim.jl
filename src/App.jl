"""
Functions to run quantum error correction simulations and merge output data.
"""
module App

# imports
using ..Model
using ..PauliTools: bsp, pack, weight
using Random: AbstractRNG, GLOBAL_RNG, MersenneTwister, from_seed, make_seed

# exports
export qec_run_once, qec_run

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
    data = Dict(
        :error_weight => weight(error),
        :logical_commutations => l_commutations,
        :success => success,
    )
    return data
end

function qec_run(
    code::StabilizerCode,
    error_model::ErrorModel,
    decoder::Decoder,
    p::Float64;
    max_runs::Union{Int, Nothing}=nothing,
    max_failures::Union{Int, Nothing}=nothing,
    random_seed=nothing
)
    @info "qec_run: starting" code=code error_model=error_model decoder=decoder p=p #=
        =# max_runs=max_runs max_failures=max_failures random_seed=random_seed
    # derived defaults
    max_runs = isnothing(max_runs) && isnothing(max_failures) ? 1 : max_runs
    wall_time_start = time_ns()

    runs_data = Dict(
        :code => label(code),
        :n_k_d => nkd(code),
        :time_steps => 1,
        :error_model => label(error_model),
        :decoder => label(decoder),
        :error_probability => p,
        :measurement_error_probability => 0.0,
        :n_run => 0,
        :n_success => 0,
        :n_fail => 0,
        :n_logical_commutations => nothing,
        :custom_totals => nothing,
        :error_weight_total => 0,
        :error_weight_pvar => 0.0,
        :logical_failure_rate => 0.0,
        :physical_error_rate => 0.0,
        :wall_time => 0.0,
    )

    rng = MersenneTwister(random_seed)
    @info "qec_run: rng=$(rng)"


    runs_data[:wall_time] = (time_ns() - wall_time_start) / 10^9

    @info "qec_run: complete: data=$(runs_data)"

    return runs_data
end

end
