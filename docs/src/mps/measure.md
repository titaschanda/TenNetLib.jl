# Measuring the MPS

```@docs
TenNetLib.entropy(psi::MPS, bond::Int)
TenNetLib.entropy(psi::MPS; kwargs...)
bond_spectrum(psi::MPS, bond::Int; by_charge = false)
bond_spectrum(psi::MPS; bonds::Union{Nothing, Vector{Int}} = nothing, by_charge = false)
measure(::Type{T}, psi::MPS, opten::ITensor) where T <: Union{ComplexF64, Float64}
measure(::Type{T}, psi::MPS, opstr::String, pos::Int) where T <: Union{ComplexF64, Float64}
measure(::Type{T}, psi::MPS, opstr::String; kwargs...) where T <: Union{ComplexF64, Float64}
measure(::Type{T}, psi::MPS, optens::Vector{ITensor}) where T <: Union{ComplexF64, Float64}
measure(::Type{T}, psi::MPS, oppairs::Vector{Pair{String, Int}}; isfermions::Bool = true) where T <: Union{ComplexF64, Float64}
```
