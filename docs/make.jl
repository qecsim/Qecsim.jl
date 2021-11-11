using Qecsim
using Documenter

DocMeta.setdocmeta!(Qecsim, :DocTestSetup, :(using Qecsim); recursive=true)

makedocs(;
    modules=[Qecsim],
    authors="David K. Tuckett, Terry Farrelly",
    repo="https://github.com/qecsim/Qecsim.jl/blob/{commit}{path}#{line}",
    sitename="Qecsim.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://qecsim.github.io/Qecsim.jl",
        assets=String[],
    ),
    pages=[
        "Overview" => "index.md",
        "API" => [
            "api/core.md",
            "api/models.md",
            "api/index.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/qecsim/Qecsim.jl",
    devbranch = "main",
)
