"""
    QecsimError(msg::String)

Indicates an internal (core or models) error.
"""
struct QecsimError <: Exception
    msg::String
end
Base.showerror(io::IO, ex::QecsimError) = print(io, "QecsimError: ", ex.msg)
