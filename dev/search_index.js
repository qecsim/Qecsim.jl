var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/#Qecsim","page":"API","title":"Qecsim","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Qecsim","category":"page"},{"location":"api/#Qecsim","page":"API","title":"Qecsim","text":"Package for simulating quantum error correction using stabilizer codes.\n\n\n\n\n\n","category":"module"},{"location":"api/","page":"API","title":"API","text":"Modules = [Qecsim]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"api/#Qecsim.doubler-Tuple{Number}","page":"API","title":"Qecsim.doubler","text":"doubler(x)\n\nReturn 2x.\n\nExamples\n\njulia> doubler(3)\n6\n\n\n\n\n\n","category":"method"},{"location":"api/#Model","page":"API","title":"Model","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"CurrentModule = Qecsim.Model","category":"page"},{"location":"api/","page":"API","title":"API","text":"Model","category":"page"},{"location":"api/#Qecsim.Model","page":"API","title":"Qecsim.Model","text":"Abstract types and methods for codes, error models and decoders.\n\n\n\n\n\n","category":"module"},{"location":"api/#Model.StabilizerCode","page":"API","title":"Model.StabilizerCode","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Model.StabilizerCode\nModel.label\nModel.logical_xs\nModel.logical_zs\nModel.logicals\nModel.nkd\nModel.stabilizers\nModel.validate","category":"page"},{"location":"api/#Qecsim.Model.StabilizerCode","page":"API","title":"Qecsim.Model.StabilizerCode","text":"Abstract supertype for stabilizer codes.\n\n\n\n\n\n","category":"type"},{"location":"api/#Qecsim.Model.label","page":"API","title":"Qecsim.Model.label","text":"label(code::StabilizerCode) -> String\n\nReturn a label suitable for use in plots and for grouping results.\n\nnote: Abstract method\nThis method should be implemented for concrete subtypes.\n\n\n\n\n\n","category":"function"},{"location":"api/#Qecsim.Model.logical_xs","page":"API","title":"Qecsim.Model.logical_xs","text":"logical_xs(code::StabilizerCode) -> BitMatrix\n\nReturn the logical X operators in binary symplectic form.\n\nEach row is an operator. The order should match that of logical_zs.\n\nnote: Abstract method\nThis method should be implemented for concrete subtypes.\n\n\n\n\n\n","category":"function"},{"location":"api/#Qecsim.Model.logical_zs","page":"API","title":"Qecsim.Model.logical_zs","text":"logical_zs(code::StabilizerCode) -> BitMatrix\n\nReturn the logical Z operators in binary symplectic form.\n\nEach row is an operator. The order should match that of logical_xs.\n\nnote: Abstract method\nThis method should be implemented for concrete subtypes.\n\n\n\n\n\n","category":"function"},{"location":"api/#Qecsim.Model.logicals","page":"API","title":"Qecsim.Model.logicals","text":"logicals(code::StabilizerCode) -> BitMatrix\n\nReturn the logical operators in binary symplectic form.\n\nEach row is an operator. X operators are stacked above Z operators in the order given by logical_xs and logical_zs.\n\n\n\n\n\n","category":"function"},{"location":"api/#Qecsim.Model.nkd","page":"API","title":"Qecsim.Model.nkd","text":"nkd(code::StabilizerCode) -> Tuple{Int, Int, Union{Int,Nothing}}\n\nReturn a descriptor in the format (n, k, d), where n is the number of physical qubits, k is the number of logical qubits, and d is the distance of the code (or nothing if unknown).\n\nnote: Abstract method\nThis method should be implemented for concrete subtypes.\n\n\n\n\n\n","category":"function"},{"location":"api/#Qecsim.Model.stabilizers","page":"API","title":"Qecsim.Model.stabilizers","text":"stabilizers(code::StabilizerCode) -> BitMatrix\n\nReturn the stabilizers in binary symplectic form.\n\nEach row is a stabilizer generator. An overcomplete set of generators can be included to simplify decoding.\n\nnote: Abstract method\nThis method should be implemented for concrete subtypes.\n\n\n\n\n\n","category":"function"},{"location":"api/#Qecsim.Model.validate","page":"API","title":"Qecsim.Model.validate","text":"validate(code::StabilizerCode)\n\nPerform sanity checks.\n\nAn AssertionError is thrown if any of the following fail:\n\nS odot S^T = 0\nS odot L^T = 0\nL odot L^T = Lambda\n\nwhere S and L are the code stabilizers and logicals, respectively, and odot and Lambda are defined in PT.bsp.\n\n\n\n\n\n","category":"function"},{"location":"api/#Model.ErrorModel","page":"API","title":"Model.ErrorModel","text":"","category":"section"},{"location":"api/#PauliTools","page":"API","title":"PauliTools","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"CurrentModule = Qecsim.PauliTools","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [PauliTools]","category":"page"},{"location":"api/#Qecsim.PauliTools","page":"API","title":"Qecsim.PauliTools","text":"Tools for Pauli strings and binary symplectic vectors / matrices.\n\n\n\n\n\n","category":"module"},{"location":"api/#Qecsim.PauliTools.bsf_to_pauli-Tuple{AbstractVector{Bool}}","page":"API","title":"Qecsim.PauliTools.bsf_to_pauli","text":"bsf_to_pauli(bsf::Union{BitVector, BitMatrix})\n\nConvert binary symplectic form to Pauli operator(s).\n\nA vector is converted to a single Pauli string. A matrix is converted row-by-row to a collection of Pauli strings.\n\nExamples\n\njulia> bsf_to_pauli(BitVector([1, 0, 0, 0, 1, 0, 0, 1, 0, 1]))\n\"XIZIY\"\n\njulia> bsf_to_pauli(BitMatrix([1 0 0 0 1 0 0 1 0 1; 0 1 0 1 0 0 0 1 1 0]))\n2-element Vector{String}:\n \"XIZIY\"\n \"IXZYI\"\n\n\n\n\n\n","category":"method"},{"location":"api/#Qecsim.PauliTools.bsp-Tuple{AbstractVecOrMat{Bool}, AbstractVecOrMat{Bool}}","page":"API","title":"Qecsim.PauliTools.bsp","text":"bsp(A::Union{BitVector, BitMatrix}, B::Union{BitVector, BitMatrix})\n\nReturn the binary symplectic product of A with B, given in binary symplectic form.\n\nThe binary symplectic product odot is defined as A odot B equiv A Lambda B bmod 2 where Lambda = leftbeginsmallmatrix 0  I  I  0 endsmallmatrixright.\n\nExamples\n\njulia> a = BitVector([1, 0, 0, 0]);  # XI\n\njulia> b = BitVector([0, 0, 1, 0]);  # ZI\n\njulia> bsp(a', b)\n1\n\njulia> stabilizers = BitMatrix(  # 5-qubit stabilizers\n       [1 0 0 1 0 0 1 1 0 0     # XZZXI\n        0 1 0 0 1 0 0 1 1 0     # IXZZX\n        1 0 1 0 0 0 0 0 1 1     # XIXZZ\n        0 1 0 1 0 1 0 0 0 1]);  # ZXIXZ\n\njulia> error = BitVector([0, 0, 1, 1, 0, 0, 1, 0, 1, 0]);  # IZXYI\n\njulia> bsp(stabilizers, error)\n4-element Vector{Int64}:\n 0\n 1\n 1\n 0\n\n\n\n\n\n","category":"method"},{"location":"api/#Qecsim.PauliTools.pauli_to_bsf-Tuple{AbstractString}","page":"API","title":"Qecsim.PauliTools.pauli_to_bsf","text":"pauli_to_bsf(pauli::Union{String, Iterable of String})\n\nConvert the Pauli operator(s) to binary symplectic form.\n\nA single Pauli string is converted to a vector. A collection of Pauli strings is converted to a matrix where each row corresponds to a Pauli.\n\nExamples\n\njulia> pauli_to_bsf(\"XIZIY\")\n10-element BitVector:\n 1\n 0\n 0\n 0\n 1\n 0\n 0\n 1\n 0\n 1\n\njulia> pauli_to_bsf([\"XIZIY\", \"IXZYI\"])\n2×10 BitMatrix:\n 1  0  0  0  1  0  0  1  0  1\n 0  1  0  1  0  0  0  1  1  0\n\n\n\n\n\n","category":"method"},{"location":"api/#Models","page":"API","title":"Models","text":"","category":"section"},{"location":"api/#BasicModels","page":"API","title":"BasicModels","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"CurrentModule = Qecsim.BasicModels","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [BasicModels]","category":"page"},{"location":"api/#Qecsim.BasicModels","page":"API","title":"Qecsim.BasicModels","text":"Basic stabilizer codes\n\n\n\n\n\n","category":"module"},{"location":"api/#Qecsim.BasicModels.BasicCode","page":"API","title":"Qecsim.BasicModels.BasicCode","text":"BasicCode <: StabilizerCode\n\nBasic code defined by its stabilizers and logical operators.\n\n\n\n\n\n","category":"type"},{"location":"api/#Qecsim.BasicModels.BasicCode-Tuple{Any, Any, Any}","page":"API","title":"Qecsim.BasicModels.BasicCode","text":"BasicCode(pauli_stabilizers, pauli_logical_xs, pauli_logical_zs; nkd=nothing,\n    label=nothing)\n\nConstruct a basic code from string representations of stabilizers and logical operators.\n\nPaulis are expressed as strings of capitalized I, X, Y, Z characters, with one character per physical qubit. Logical X and Z operators are in matching order, with one of each for each logical qubit. Optional nkd defaults to n and k evaluated and d nothing. Optional label defaults to \"Basic [n, k, d]\".\n\nExamples\n\njulia> using Qecsim.BasicModels\n\njulia> code = BasicCode([\"ZZI\", \"IZZ\"], [\"XXX\"], [\"IIZ\"]);  # 3-qubit repetition\n\njulia> validate(code)  # no error indicates operators satisfy commutation relations\n\njulia> nkd(code)  # default nkd\n(3, 1, nothing)\n\njulia> label(code)  # default label\n\"Basic [3,1,nothing]\"\n\n\n\n\n\n","category":"method"},{"location":"api/#Qecsim.BasicModels.FiveQubitCode-Tuple{}","page":"API","title":"Qecsim.BasicModels.FiveQubitCode","text":"FiveQubitCode()\n\nConstruct 5-qubit [5,1,3] code as a BasicCode.\n\n\n\n\n\n","category":"method"},{"location":"api/#Qecsim.BasicModels.SteaneCode-Tuple{}","page":"API","title":"Qecsim.BasicModels.SteaneCode","text":"SteaneCode()\n\nConstruct Steane [7,1,3] code as a BasicCode.\n\n\n\n\n\n","category":"method"},{"location":"api_index/#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"api_index/","page":"Index","title":"Index","text":"","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"CurrentModule = Qecsim","category":"page"},{"location":"#Qecsim.jl-Quantum-Error-Correction-Simulator","page":"Overview","title":"Qecsim.jl - Quantum Error Correction Simulator","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"Qecsim.jl is a Julia package for simulating quantum error correction using stabilizer codes.","category":"page"},{"location":"#Contents","page":"Overview","title":"Contents","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"Pages = [\"api.md\", \"api_index.md\"]","category":"page"}]
}
