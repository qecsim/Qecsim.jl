var documenterSearchIndex = {"docs":
[{"location":"api/core/#Core","page":"Core","title":"Core","text":"","category":"section"},{"location":"api/core/#Qecsim","page":"Core","title":"Qecsim","text":"","category":"section"},{"location":"api/core/","page":"Core","title":"Core","text":"Qecsim\nQecsim.QecsimError","category":"page"},{"location":"api/core/#Qecsim","page":"Core","title":"Qecsim","text":"Package for simulating quantum error correction using stabilizer codes.\n\n\n\n\n\n","category":"module"},{"location":"api/core/#Qecsim.QecsimError","page":"Core","title":"Qecsim.QecsimError","text":"QecsimError <: Exception\n\nQecsimError(msg)\n\nConstruct an exception indicating an internal (core or models) error.\n\n\n\n\n\n","category":"type"},{"location":"api/core/#Model","page":"Core","title":"Model","text":"","category":"section"},{"location":"api/core/","page":"Core","title":"Core","text":"CurrentModule = Qecsim.Model","category":"page"},{"location":"api/core/","page":"Core","title":"Core","text":"Model","category":"page"},{"location":"api/core/#Qecsim.Model","page":"Core","title":"Qecsim.Model","text":"Abstract types and methods for codes, error models and decoders.\n\n\n\n\n\n","category":"module"},{"location":"api/core/#Model.AbstractModel","page":"Core","title":"Model.AbstractModel","text":"","category":"section"},{"location":"api/core/","page":"Core","title":"Core","text":"Model.AbstractModel\nModel.label","category":"page"},{"location":"api/core/#Qecsim.Model.AbstractModel","page":"Core","title":"Qecsim.Model.AbstractModel","text":"Abstract supertype for models.\n\n\n\n\n\n","category":"type"},{"location":"api/core/#Qecsim.Model.label","page":"Core","title":"Qecsim.Model.label","text":"label(code::AbstractModel) -> String\n\nReturn a label suitable for use in plots and for grouping results.\n\nnote: Abstract method\nThis method should be implemented for concrete subtypes of AbstractModel.\n\n\n\n\n\n","category":"function"},{"location":"api/core/#Model.StabilizerCode","page":"Core","title":"Model.StabilizerCode","text":"","category":"section"},{"location":"api/core/","page":"Core","title":"Core","text":"Model.StabilizerCode\nModel.logical_xs\nModel.logical_zs\nModel.logicals\nModel.nkd\nModel.stabilizers\nModel.validate","category":"page"},{"location":"api/core/#Qecsim.Model.StabilizerCode","page":"Core","title":"Qecsim.Model.StabilizerCode","text":"StabilizerCode <: AbstractModel\n\nAbstract supertype for stabilizer codes.\n\n\n\n\n\n","category":"type"},{"location":"api/core/#Qecsim.Model.logical_xs","page":"Core","title":"Qecsim.Model.logical_xs","text":"logical_xs(code::StabilizerCode) -> BitMatrix\n\nReturn the logical X operators in binary symplectic form.\n\nEach row is an operator. The order should match that of logical_zs.\n\nnote: Abstract method\nThis method should be implemented for concrete subtypes of StabilizerCode.\n\n\n\n\n\n","category":"function"},{"location":"api/core/#Qecsim.Model.logical_zs","page":"Core","title":"Qecsim.Model.logical_zs","text":"logical_zs(code::StabilizerCode) -> BitMatrix\n\nReturn the logical Z operators in binary symplectic form.\n\nEach row is an operator. The order should match that of logical_xs.\n\nnote: Abstract method\nThis method should be implemented for concrete subtypes of StabilizerCode.\n\n\n\n\n\n","category":"function"},{"location":"api/core/#Qecsim.Model.logicals","page":"Core","title":"Qecsim.Model.logicals","text":"logicals(code::StabilizerCode) -> BitMatrix\n\nReturn the logical operators in binary symplectic form.\n\nEach row is an operator. X operators are stacked above Z operators in the order given by logical_xs and logical_zs.\n\n\n\n\n\n","category":"function"},{"location":"api/core/#Qecsim.Model.nkd","page":"Core","title":"Qecsim.Model.nkd","text":"nkd(code::StabilizerCode) -> Tuple{Int, Int, Union{Int, Missing}}\n\nReturn a descriptor in the format (n, k, d), where n is the number of physical qubits, k is the number of logical qubits, and d is the distance of the code (or missing if unknown).\n\nnote: Abstract method\nThis method should be implemented for concrete subtypes of StabilizerCode.\n\n\n\n\n\n","category":"function"},{"location":"api/core/#Qecsim.Model.stabilizers","page":"Core","title":"Qecsim.Model.stabilizers","text":"stabilizers(code::StabilizerCode) -> BitMatrix\n\nReturn the stabilizers in binary symplectic form.\n\nEach row is a stabilizer generator. An overcomplete set of generators can be included to simplify decoding.\n\nnote: Abstract method\nThis method should be implemented for concrete subtypes of StabilizerCode.\n\n\n\n\n\n","category":"function"},{"location":"api/core/#Qecsim.Model.validate","page":"Core","title":"Qecsim.Model.validate","text":"validate(code::StabilizerCode)\n\nPerform sanity checks.\n\nIf any of the following fail then a QecsimError is thrown:\n\nS odot S^T = 0\nS odot L^T = 0\nL odot L^T = Lambda\n\nwhere S and L are the code stabilizers and logicals, respectively, and odot and Lambda are defined in bsp.\n\n\n\n\n\n","category":"function"},{"location":"api/core/#Model.ErrorModel","page":"Core","title":"Model.ErrorModel","text":"","category":"section"},{"location":"api/core/","page":"Core","title":"Core","text":"Model.ErrorModel\nModel.generate\nModel.probability_distribution","category":"page"},{"location":"api/core/#Qecsim.Model.ErrorModel","page":"Core","title":"Qecsim.Model.ErrorModel","text":"ErrorModel <: AbstractModel\n\nAbstract supertype for error models.\n\n\n\n\n\n","category":"type"},{"location":"api/core/#Qecsim.Model.generate","page":"Core","title":"Qecsim.Model.generate","text":"generate(error_model::ErrorModel, code::StabilizerCode, p::Float64,\n         [rng::AbstractRNG=GLOBAL_RNG]) -> BitVector\n\nGenerate a new error in binary symplectic form according to the error_model and code, where p is typically the probability of an error on a single qubit.\n\nnote: Abstract method\nThis method should be implemented for concrete subtypes of ErrorModel.\n\n\n\n\n\n","category":"function"},{"location":"api/core/#Qecsim.Model.probability_distribution","page":"Core","title":"Qecsim.Model.probability_distribution","text":"probability_distribution(error_model::ErrorModel, p::Float64) -> NTuple{4, Real}\n\nReturn the single-qubit probability distribution amongst Pauli I, X, Y and Z, where p is the overall probability of an error on a single qubit.\n\nnote: Abstract method [optional]\nThis method is not invoked by any core modules. Since it is often useful for decoders, it is provided as a template and concrete subtypes of ErrorModel are encouraged to implement it when appropriate, particularly for IID error models.\n\n\n\n\n\n","category":"function"},{"location":"api/core/#Model.Decoder","page":"Core","title":"Model.Decoder","text":"","category":"section"},{"location":"api/core/","page":"Core","title":"Core","text":"Model.Decoder\nModel.DecodeResult\nModel.decode","category":"page"},{"location":"api/core/#Qecsim.Model.Decoder","page":"Core","title":"Qecsim.Model.Decoder","text":"Decoder <: AbstractModel\n\nAbstract supertype for decoders.\n\n\n\n\n\n","category":"type"},{"location":"api/core/#Qecsim.Model.DecodeResult","page":"Core","title":"Qecsim.Model.DecodeResult","text":"DecodeResult(recovery::AbstractVector{Bool})\n\nConstruct a decoding result including the recovery operation.\n\nwarning: Warning\nIn the future, this type will be extended to include success flags and more.\n\n\n\n\n\n","category":"type"},{"location":"api/core/#Qecsim.Model.decode","page":"Core","title":"Qecsim.Model.decode","text":"decode(decoder::Decoder, code::StabilizerCode, syndrome::AbstractVector{Bool};\n       kwargs...) -> DecodeResult\n\nResolve a recovery operation for the given code and syndrome, or evaluate the success of decoding, as encapsulated in the decode result.\n\nThe syndrome has length equal to the number of code stabilizers, and element values of 0 or 1 (false or true) indicate whether the corresponding stabilizer does or does not commute with the error, respectively.\n\nKeyword parameters kwargs may be provided by the client with context values such as error_model, error_probability and error. Most implementations will ignore such parameters; however, if they are used, implementations should declare them explicitly and treat them as optional.\n\nnote: Abstract method\nThis method should be implemented for concrete subtypes of Decoder.\n\n\n\n\n\n","category":"function"},{"location":"api/core/#PauliTools","page":"Core","title":"PauliTools","text":"","category":"section"},{"location":"api/core/","page":"Core","title":"Core","text":"CurrentModule = Qecsim.PauliTools","category":"page"},{"location":"api/core/","page":"Core","title":"Core","text":"PauliTools\nPauliTools.bsp\nPauliTools.to_bsf\nPauliTools.to_pauli\nPauliTools.weight","category":"page"},{"location":"api/core/#Qecsim.PauliTools","page":"Core","title":"Qecsim.PauliTools","text":"Tools for Pauli strings and binary symplectic vectors / matrices.\n\n\n\n\n\n","category":"module"},{"location":"api/core/#Qecsim.PauliTools.bsp","page":"Core","title":"Qecsim.PauliTools.bsp","text":"bsp(A::AbstractVecOrMat{Bool}, B::AbstractVecOrMat{Bool})\n    -> Union{Bool, AbstractVecOrMat{Bool}}\n\nReturn the binary symplectic product of A with B, given in binary symplectic form.\n\nThe binary symplectic product odot is defined as A odot B equiv A Lambda B bmod 2 where Lambda = leftbeginsmallmatrix 0  I  I  0 endsmallmatrixright.\n\nExamples\n\njulia> a = BitVector([1, 0, 0, 0]);  # XI\n\njulia> b = BitVector([0, 0, 1, 0]);  # ZI\n\njulia> bsp(a', b)\ntrue\n\njulia> stabilizers = BitMatrix(  # 5-qubit stabilizers\n       [1 0 0 1 0 0 1 1 0 0     # XZZXI\n        0 1 0 0 1 0 0 1 1 0     # IXZZX\n        1 0 1 0 0 0 0 0 1 1     # XIXZZ\n        0 1 0 1 0 1 0 0 0 1]);  # ZXIXZ\n\njulia> error = BitVector([0, 0, 1, 1, 0, 0, 1, 0, 1, 0]);  # IZXYI\n\njulia> bsp(stabilizers, error)\n4-element BitVector:\n 0\n 1\n 1\n 0\n\n\n\n\n\n","category":"function"},{"location":"api/core/#Qecsim.PauliTools.to_bsf","page":"Core","title":"Qecsim.PauliTools.to_bsf","text":"to_bsf(pauli::Union{AbstractString, AbstractVector{<:AbstractString}})\n\nConvert the Pauli string operator(s) to binary symplectic form.\n\nA single Pauli string is converted to a vector. A vector of Pauli strings is converted to a matrix where each row corresponds to a Pauli.\n\nExamples\n\njulia> to_bsf(\"XIZIY\")\n10-element BitVector:\n 1\n 0\n 0\n 0\n 1\n 0\n 0\n 1\n 0\n 1\n\njulia> to_bsf([\"XIZIY\", \"IXZYI\"])\n2×10 BitMatrix:\n 1  0  0  0  1  0  0  1  0  1\n 0  1  0  1  0  0  0  1  1  0\n\n\n\n\n\n","category":"function"},{"location":"api/core/#Qecsim.PauliTools.to_pauli","page":"Core","title":"Qecsim.PauliTools.to_pauli","text":"to_pauli(bsf::AbstractVecOrMat{Bool})\n\nConvert the binary symplectic form to Pauli string operator(s).\n\nA vector is converted to a single Pauli string. A matrix is converted row-by-row to a collection of Pauli strings.\n\nExamples\n\njulia> to_pauli(BitVector([1, 0, 0, 0, 1, 0, 0, 1, 0, 1]))\n\"XIZIY\"\n\njulia> to_pauli(BitMatrix([1 0 0 0 1 0 0 1 0 1; 0 1 0 1 0 0 0 1 1 0]))\n2-element Vector{String}:\n \"XIZIY\"\n \"IXZYI\"\n\n\n\n\n\n","category":"function"},{"location":"api/core/#Qecsim.PauliTools.weight","page":"Core","title":"Qecsim.PauliTools.weight","text":"weight(bsf::AbstractVecOrMat{Bool})\n\nReturn the weight of the binary symplectic form.\n\nExamples\n\njulia> weight(BitVector([1, 0, 0, 0, 1, 0, 0, 1, 0, 1]))\n3\n\njulia> weight(BitMatrix([1 0 0 0 1 0 0 1 0 1; 1 1 1 1 1 0 0 0 0 0]))\n8\n\n\n\n\n\nweight(pauli::Union{AbstractString, AbstractVector{<:AbstractString}})\n\nReturn the weight of the Pauli string operator(s).\n\nExamples\n\njulia> weight(\"XIZIY\")\n3\n\njulia> weight([\"XIZIY\", \"XXXXX\"])\n8\n\n\n\n\n\n","category":"function"},{"location":"api/#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"api/","page":"Index","title":"Index","text":"","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"CurrentModule = Qecsim","category":"page"},{"location":"#Qecsim.jl-Quantum-Error-Correction-Simulator","page":"Overview","title":"Qecsim.jl - Quantum Error Correction Simulator","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"Qecsim.jl is a Julia package for simulating quantum error correction using stabilizer codes.","category":"page"},{"location":"#Contents","page":"Overview","title":"Contents","text":"","category":"section"},{"location":"#API","page":"Overview","title":"API","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"Pages = [\"api/core.md\", \"api/models.md\", \"api/index.md\"]","category":"page"},{"location":"api/models/#Models","page":"Models","title":"Models","text":"","category":"section"},{"location":"api/models/#BasicModels","page":"Models","title":"BasicModels","text":"","category":"section"},{"location":"api/models/","page":"Models","title":"Models","text":"CurrentModule = Qecsim.BasicModels","category":"page"},{"location":"api/models/","page":"Models","title":"Models","text":"BasicModels\nBasicModels.BasicCode\nBasicModels.FiveQubitCode\nBasicModels.SteaneCode","category":"page"},{"location":"api/models/#Qecsim.BasicModels","page":"Models","title":"Qecsim.BasicModels","text":"Basic stabilizer codes\n\n\n\n\n\n","category":"module"},{"location":"api/models/#Qecsim.BasicModels.BasicCode","page":"Models","title":"Qecsim.BasicModels.BasicCode","text":"BasicCode <: StabilizerCode\n\nBasicCode(pauli_stabilizers, pauli_logical_xs, pauli_logical_zs;\n          nkd=nothing, label=nothing)\n\nConstruct a basic code from string representations of stabilizers and logical operators.\n\nPaulis are expressed as strings of capitalized I, X, Y, Z characters, with one character per physical qubit. Logical X and Z operators are in matching order, with one of each for each logical qubit. Optional nkd defaults to n and k evaluated and d missing. Optional label defaults to \"Basic [n,k,d]\".\n\nExamples\n\njulia> using Qecsim.BasicModels\n\njulia> code = BasicCode([\"ZZI\", \"IZZ\"], [\"XXX\"], [\"IIZ\"]);  # 3-qubit repetition\n\njulia> validate(code)  # no error indicates operators satisfy commutation relations\n\njulia> nkd(code)  # default nkd\n(3, 1, missing)\n\njulia> label(code)  # default label\n\"Basic [3,1,missing]\"\n\n\n\n\n\n","category":"type"},{"location":"api/models/#Qecsim.BasicModels.FiveQubitCode","page":"Models","title":"Qecsim.BasicModels.FiveQubitCode","text":"FiveQubitCode() -> BasicCode\n\nConstruct 5-qubit [5,1,3] code as a BasicCode.\n\n\n\n\n\n","category":"function"},{"location":"api/models/#Qecsim.BasicModels.SteaneCode","page":"Models","title":"Qecsim.BasicModels.SteaneCode","text":"SteaneCode() -> BasicCode\n\nConstruct Steane [7,1,3] code as a BasicCode.\n\n\n\n\n\n","category":"function"},{"location":"api/models/#GenericModels","page":"Models","title":"GenericModels","text":"","category":"section"},{"location":"api/models/","page":"Models","title":"Models","text":"CurrentModule = Qecsim.GenericModels","category":"page"},{"location":"api/models/","page":"Models","title":"Models","text":"GenericModels\nGenericModels.SimpleErrorModel\nModel.generate(::SimpleErrorModel, ::StabilizerCode, ::Float64, ::AbstractRNG)\nGenericModels.BitFlipErrorModel\nGenericModels.BitPhaseFlipErrorModel\nGenericModels.DepolarizingErrorModel\nGenericModels.PhaseFlipErrorModel\nGenericModels.NaiveDecoder","category":"page"},{"location":"api/models/#Qecsim.GenericModels","page":"Models","title":"Qecsim.GenericModels","text":"Generic error models and decoders compatible with any stabilizer codes.\n\n\n\n\n\n","category":"module"},{"location":"api/models/#Qecsim.GenericModels.SimpleErrorModel","page":"Models","title":"Qecsim.GenericModels.SimpleErrorModel","text":"SimpleErrorModel <: ErrorModel\n\nAbstract supertype for simple IID error models that generate errors based on the number of qubits and a probability distribution.\n\n\n\n\n\n","category":"type"},{"location":"api/models/#Qecsim.Model.generate-Tuple{Qecsim.GenericModels.SimpleErrorModel, StabilizerCode, Float64, Random.AbstractRNG}","page":"Models","title":"Qecsim.Model.generate","text":"generate(error_model::SimpleErrorModel, code::StabilizerCode, p::Float64,\n         [rng::AbstractRNG=GLOBAL_RNG]) -> BitVector\n\nGenerate a new IID error based on probability_distribution. See also Model.generate.\n\nnote: Note\nThe method probability_distribution should be implemented for concrete subtypes of SimpleErrorModel.\n\n\n\n\n\n","category":"method"},{"location":"api/models/#Qecsim.GenericModels.BitFlipErrorModel","page":"Models","title":"Qecsim.GenericModels.BitFlipErrorModel","text":"BitFlipErrorModel <: SimpleErrorModel\n\nIID error model with probability vector: (p_I p_X p_Y p_Z) = (1-p p 0 0), where p is the probability of an error on a single-qubit.\n\n\n\n\n\n","category":"type"},{"location":"api/models/#Qecsim.GenericModels.BitPhaseFlipErrorModel","page":"Models","title":"Qecsim.GenericModels.BitPhaseFlipErrorModel","text":"BitPhaseFlipErrorModel <: SimpleErrorModel\n\nIID error model with probability vector: (p_I p_X p_Y p_Z) = (1-p 0 p 0), where p is the probability of an error on a single-qubit.\n\n\n\n\n\n","category":"type"},{"location":"api/models/#Qecsim.GenericModels.DepolarizingErrorModel","page":"Models","title":"Qecsim.GenericModels.DepolarizingErrorModel","text":"DepolarizingErrorModel <: SimpleErrorModel\n\nIID error model with probability vector: (p_I p_X p_Y p_Z) = (1-p p3 p3 p3), where p is the probability of an error on a single-qubit.\n\n\n\n\n\n","category":"type"},{"location":"api/models/#Qecsim.GenericModels.PhaseFlipErrorModel","page":"Models","title":"Qecsim.GenericModels.PhaseFlipErrorModel","text":"PhaseFlipErrorModel <: SimpleErrorModel\n\nIID error model with probability vector: (p_I p_X p_Y p_Z) = (1-p 0 0 p), where p is the probability of an error on a single-qubit.\n\n\n\n\n\n","category":"type"},{"location":"api/models/#Qecsim.GenericModels.NaiveDecoder","page":"Models","title":"Qecsim.GenericModels.NaiveDecoder","text":"NaiveDecoder <: Decoder\n\nNaiveDecoder([max_qubits=10])\n\nConstruct a naive decoder that iterates through all possible errors, in ascending weight, and resolves to the first error that matches the syndrome.\n\nnote: Note\nThis decoder is slow for even moderate numbers of qubits. By default, it is restricted to codes with a maximum of 10 qubits.\n\n\n\n\n\n","category":"type"}]
}
