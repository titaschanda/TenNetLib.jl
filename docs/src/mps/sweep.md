# Sweeping through the MPS

At a lower level of abstraction, TeNLib.jl allows to control each fullsweep
(left-to-right and right-to-left) manually to update `StateEnvs`.

Skip this part if you want to avoid lower-level abstraction.

## `SweepData`

TeNLib.jl defines a struct, called  `SweepData`, to store essential data after each fullsweep.

```@docs
SweepData
Base.copy(swdata::SweepData)
```

## Perform a fullsweep

```@docs
fullsweep!(sysenv::StateEnvs, solver, nsite::Int, swdata::SweepData; kwargs...)
```

## Perform a dynamical fullsweep

TeNLib.jl defines the following function to dynamically decide whether to perform single- or
two-site update at each bond, depending on the entropy growth at the previous halfsweep.

```@docs
dynamic_fullsweep!(sysenv::StateEnvs, solver, swdata::SweepData; kwargs...)
```

## Global Subspace Expansion

Following [Phys. Rev. B **102**, 094315 (2020)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.094315), a Global Subspace Expansion can be performed using Krylov subspace if the
environments are created by a single `MPO`.

Apart from TDVP, Global Subspace Expansion is also very useful for DMRG to get rid of nasty
local minimas.

```@docs
krylov_extend!(psi::MPS, H::MPO; kwargs...)
krylov_extend!(sysenv::StateEnvs{ProjMPO}; kwargs...)
```

