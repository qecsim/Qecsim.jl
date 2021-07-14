using Test
using Qecsim: QecsimError

# taken from https://github.com/JuliaLang/julia/blob/v1.6.1/test/errorshow.jl
macro except_str(expr, err_type)
    return quote
        let err = nothing
            try
                $(esc(expr))
            catch err
            end
            err === nothing && error("expected failure, but no exception thrown")
            @test typeof(err) === $(esc(err_type))
            buf = IOBuffer()
            showerror(buf, err)
            String(take!(buf))
        end
    end
end

@testset "QecsimError" begin
    @test_throws QecsimError throw(QecsimError("something went wrong!"))
    err_str = @except_str throw(QecsimError("something went wrong!")) QecsimError
    @test err_str == "QecsimError: something went wrong!"
end
