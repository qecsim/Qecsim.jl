"""
Functions to run quantum error correction simulations and merge output data.
"""
module App

# imports
using ..Model
using ..PauliTools: bsp, pack, weight
using Random: AbstractRNG, GLOBAL_RNG

# exports
export qec_run_once

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

end
