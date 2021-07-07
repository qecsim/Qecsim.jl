# Models

## BasicModels
```@meta
CurrentModule = Qecsim.BasicModels
```
```@docs
BasicModels
BasicModels.BasicCode
BasicModels.FiveQubitCode
BasicModels.SteaneCode
```

## GenericModels
```@meta
CurrentModule = Qecsim.GenericModels
```
```@docs
GenericModels
```
### GenericModels: Error Models
```@docs
GenericModels.SimpleErrorModel
Model.generate(::SimpleErrorModel, ::StabilizerCode, ::Float64, ::AbstractRNG)
GenericModels.BitFlipErrorModel
GenericModels.BitPhaseFlipErrorModel
GenericModels.DepolarizingErrorModel
GenericModels.PhaseFlipErrorModel
```
### GenericModels: Decoders
```@docs
GenericModels.NaiveDecoder
```
