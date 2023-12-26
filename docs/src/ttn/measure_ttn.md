# Measuring the TTN

```@docs
measure(::Type{T}, psi::TTN, opten::ITensor) where T <: Union{ComplexF64, Float64}
measure(::Type{T}, psi::TTN, opstr::String, pos::Int) where T <: Union{ComplexF64, Float64}
measure(::Type{T}, psi::TTN, opstr::String; kwargs...) where T <: Union{ComplexF64, Float64}
measure(::Type{T}, psi::TTN, optens::Vector{ITensor}) where T <: Union{ComplexF64, Float64}
measure(::Type{T}, psi::TTN, oppairs::Vector{Pair{String, Int}}; isfermions::Bool = true) where T <: Union{ComplexF64, Float64}
```