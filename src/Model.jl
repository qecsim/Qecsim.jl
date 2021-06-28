module Model

using Qecsim:PauliTools as PT

"""
Supertype for stabilizer codes.
"""
abstract type StabilizerCode end

"""
    stabilizers(code::StabilizerCode) -> BitMatrix

Return `stabilizers` in binary symplectic form. Each row is a stabilizer generator. An
overcomplete set of generators can be included to simplify decoding.
"""
function stabilizers end
function logical_xs end
function logical_zs end
function logicals(code::StabilizerCode)
    vcat(logical_xs(code), logical_zs(code))
end
function nkd end
function label end
function validate(code::StabilizerCode)
    s, l = stabilizers(code), logicals(code)
    @assert all(PT.bsp(s, transpose(s)) .== 0)
    @assert all(PT.bsp(s, transpose(l)) .== 0)
    # twisted identity with same size as logicals
    # TODO: complete
end

end
