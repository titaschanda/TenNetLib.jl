# Sweeping through the TTN

At a lower level of abstraction, TenNetLib.jl allows to control each fullsweep
 manually to update `StateEnvsTTN`.

Skip this part if you want to avoid lower-level abstraction.

## `SweepDataTTN`

TenNetLib.jl defines a struct, called  `SweepDataTTN`, to store essential data after each fullsweep.

```@docs
SweepDataTTN
Base.copy(swdata::SweepDataTTN)
```
## `sweeppath`

```@docs
default_sweeppath
```

## Perform a fullsweep

```@docs
fullsweep!(sysenv::StateEnvsTTN, sweeppath::Vector{Int2}, solver, swdata::SweepDataTTN; kwargs...)
```
