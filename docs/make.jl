using Qecsim
using Documenter

DocMeta.setdocmeta!(Qecsim, :DocTestSetup, :(using Qecsim); recursive=true)

makedocs(;
    modules=[Qecsim],
    authors="David K. Tuckett, Terry Farrelly",
    repo="https://github.com/dkt29/Qecsim.jl/blob/{commit}{path}#{line}",
    sitename="Qecsim.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dkt29.github.io/Qecsim.jl",
        assets=String[],
    ),
    pages=[
        "Overview" => "index.md",
        "API" => "api.md",
        "Index" => "api_index.md"
    ],
)

deploydocs(;
    repo="github.com/dkt29/Qecsim.jl",
    devbranch = "main",
)
