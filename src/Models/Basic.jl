module Basic

import ...Model

using ...Model:StabilizerCode

struct BasicCode <: StabilizerCode
    label::String
    function BasicCode(label::Union{String, Nothing}=nothing)
        # TODO: Should label be a keyword argument with nkd?
        label = isnothing(label) ? "BasicCode" : label
        new(label)
    end
end
function Model.label(code::BasicCode)
    code.label
end

end
