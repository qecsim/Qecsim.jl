```@meta
CurrentModule = Qecsim
```

# Qecsim.jl - Quantum Error Correction Simulator

## Introduction

[Qecsim.jl](https://github.com/qecsim/Qecsim.jl) is a Julia package for
simulating quantum error correction using stabilizer codes.

**NOTE**: _Qecsim.jl_ is a ground-up rewrite of the Python package
[qecsim](https://github.com/qecsim/qecsim).

## Installation

_Qecsim.jl_ is installed, like any other registered Julia package, using the Julia package
manager [Pkg](https://pkgdocs.julialang.org/):

```julia
pkg> add Qecsim  # Press ']' to enter the Pkg REPL mode.
```
or
```julia
julia> using Pkg; Pkg.add("Qecsim")
```

## Examples

### Simulation run

```@repl
using Qecsim, Qecsim.BasicModels, Qecsim.GenericModels
data = qec_run(FiveQubitCode(), BitFlipErrorModel(), NaiveDecoder(), 0.1; max_runs=100);
data
```

### Simulation plot

**NOTE**: This example assumes that [Plots](https://docs.juliaplots.org) is installed.

```@repl
using Qecsim, Qecsim.BasicModels, Qecsim.GenericModels, Logging, Plots
with_logger(NullLogger()) do
    error_probabilities = 0.0:0.02:0.5
    codes = [FiveQubitCode(), SteaneCode()]
    error_model = DepolarizingErrorModel()
    decoder = NaiveDecoder()
    max_runs = 5000
    failure_rates = [[] for _ in codes]
    for p in error_probabilities, (code, f) in zip(codes, failure_rates)
        data = qec_run(code, error_model, decoder, p; max_runs=max_runs)
        push!(f, data[:logical_failure_rate])
    end
    labels = reshape([label(c) for c in codes], 1, :)
    plot(
        error_probabilities, failure_rates; label=labels, legend=:right,
        xlabel="Error probability", ylabel="Logical failure rate"
    )
end;
savefig("examples-plot.png"); # hide
```
![examples plot](examples-plot.png)

## Citing

Please cite _Qecsim.jl_ if you use it in your research. It was first introduced in the
following thesis:

* D. K. Tuckett,
  _Tailoring surface codes: Improvements in quantum error correction with biased noise_,
  [Ph.D. thesis](https://doi.org/10.25910/x8xw-9077),
  University of Sydney (2020), (qecsim: <https://github.com/qecsim/Qecsim.jl>).

A suitable BibTeX entry is:

    @phdthesis{qecsim,
        author = {Tuckett, David Kingsley},
        title = {Tailoring surface codes: Improvements in quantum error correction with biased noise},
        school = {University of Sydney},
        doi = {10.25910/x8xw-9077},
        year = {2020},
        note = {(qecsim: \url{https://github.com/qecsim/Qecsim.jl})}
    }

## License

_Qecsim.jl_ is released under the BSD 3-Clause license, see
[LICENSE](https://github.com/qecsim/Qecsim.jl/blob/main/LICENSE).

## Links

* Source code: <https://github.com/qecsim/Qecsim.jl>
* Documentation: <https://qecsim.github.io/Qecsim.jl>
* Contact: [qecsim@gmail.com](mailto:qecsim@gmail.com)
