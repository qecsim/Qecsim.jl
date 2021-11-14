# Qecsim.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://qecsim.github.io/Qecsim.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://qecsim.github.io/Qecsim.jl/dev/)
[![Build Status](https://github.com/qecsim/Qecsim.jl/workflows/CI/badge.svg)](https://github.com/qecsim/Qecsim.jl/actions)
[![Coverage](https://codecov.io/gh/qecsim/Qecsim.jl/branch/main/graph/badge.svg?token=nzuAO7xE6r)](https://codecov.io/gh/qecsim/Qecsim.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

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

See _Qecsim.jl_ [documentation](https://qecsim.github.io/Qecsim.jl) for more examples.

### Simulation run

```julia
julia> using Qecsim, Qecsim.BasicModels, Qecsim.GenericModels

julia> data = qec_run(FiveQubitCode(), BitFlipErrorModel(), NaiveDecoder(), 0.1; max_runs=100);

julia> data
Dict{Symbol, Any} with 17 entries:
  :error_weight_pvar             => 0.4451
  :time_steps                    => 1
  :n_logical_commutations        => [8, 3]
  :error_weight_total            => 43
  :wall_time                     => 0.00264957
  :n_k_d                         => (5, 1, 3)
  :error_model                   => "Bit-flip"
  :physical_error_rate           => 0.086
  :measurement_error_probability => 0.0
  :error_probability             => 0.1
  :n_success                     => 92
  :logical_failure_rate          => 0.08
  :custom_totals                 => nothing
  :code                          => "5-qubit"
  :decoder                       => "Naive"
  :n_fail                        => 8
  :n_run                         => 100
```

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
