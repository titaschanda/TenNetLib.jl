# Optimizing TTN

Functions to optimize TTN.

## `OptimizeParamsTTN`

TenNetLib.jl defines a struct, called  `OptimizeParamsTTN`, to control TTN optimizations.

```@docs
OptimizeParamsTTN
OptimizeParamsTTN(;maxdim::Vector{Int}, nsweeps::Vector{Int}, cutoff::Union{Vector{Float64}, Float64} = _Float64_Threshold, noise::Union{Vector{Float64}, Float64, Int} = 0.0, noisedecay::Union{Vector{Float64}, Float64, Int} = 1.0, disable_noise_after::Union{Vector{Int}, Int} = typemax(Int))
Base.copy(params::OptimizeParamsTTN)
```

## A lower level optimization function

Following function modifies `StateEnvsTTN` in-place. Skip this function if you want to
avoid lower-level abstraction.

```@docs
optimize!(sysenv::StateEnvsTTN, params::OptimizeParamsTTN,sweeppath::Vector{Int2}; kwargs...)
```

## Higher level optimzation functions

```@docs
optimize(psi0::TTN, H::CouplingModel, params::OptimizeParamsTTN, sweeppath::Vector{Int2}; kwargs...)
```


