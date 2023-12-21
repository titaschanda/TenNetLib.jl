# Measuring the MPS

```@docs
TeNLib.entropy(psi::MPS, bond::Int)
TeNLib.entropy(psi::MPS; kwargs...)
expectC(psi::MPS, opten::ITensor)
expectC(psi::MPS, opstr::String, pos::Int)
expectC(psi::MPS, opstr::String; kwargs...)
expectC(psi::MPS, optens::Vector{ITensor})
expectC(psi::MPS, oppairs::Vector{Pair{String, Int}}; isfermions::Bool = true)
```