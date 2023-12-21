# Performing DMRG

Functions to perfom (and control) DMRG runs.

## `DMRGParams`

TeNLib.jl defines a struct, called  `DMRGParams`, to control DMRG simulations.

```@docs
DMRGParams
DMRGParams(;maxdim::Vector{Int}, nsweeps::Vector{Int}, cutoff::Union{Vector{Float64}, Float64} = _Float64_Threshold, noise::Union{Vector{Float64}, Float64, Int} = 0.0, noisedecay::Union{Vector{Float64}, Float64, Int} = 1.0, disable_noise_after::Union{Vector{Int}, Int} = typemax(Int))
Base.copy(params::DMRGParams)
```

## A lower level DMRG function

Following function modifies `StateEnvs` in-place. Skip this function if you want to avoid lower-level abstraction.

```@docs
dmrg!(sysenv::StateEnvs, params::DMRGParams, nsite::Int; kwargs...)
```

## Higher level DMRG functions

Standard two- and single-site DMRG functions. Single-site DMRG can increasing the bond-dimension
if `noise > Float64_threshold()`.
```@docs
dmrg2(psi0::MPS, H::T, params::DMRGParams; kwargs...) where T <: Union{MPO, Vector{MPO}, CouplingModel}
dmrg1(psi0::MPS, H::T, params::DMRGParams; kwargs...) where T <: Union{MPO, Vector{MPO}, CouplingModel}
```

