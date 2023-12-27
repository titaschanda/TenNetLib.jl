# `StateEnvs`: A container to store the MPS and its environments

At the lowest-level of abstraction, TenNetLib.jl defines `StateEnvs` to hold an MPS and its
environments to be modified in place.

Skip this part if you want to avoid lower-level abstraction.

```@docs
StateEnvs
getpsi(sysenv::StateEnvs)
getenv(sysenv::StateEnvs)
StateEnvs(psi::MPS, H::MPO)
StateEnvs(psi::MPS, Hs::Vector{MPO})
StateEnvs(psi::MPS, H::MPO, Ms::Vector{MPS}; weight::Float64)
StateEnvs(psi::MPS, Hs::Vector{MPO}, Ms::Vector{MPS}; weight::Float64)
StateEnvs(psi::MPS, H::CouplingModel)
StateEnvs(psi::MPS, H::CouplingModel, Ms::Vector{MPS}; weight::Float64)
updateH!(sysenv::StateEnvs{ProjMPO}, H::MPO; recalcEnv::Bool = true)
updateH!(sysenv::StateEnvs{ProjMPOSum2}, Hs::Vector{MPO}; recalcEnv::Bool = true)
updateH!(sysenv::StateEnvs{ProjMPO_MPS2}, H::MPO, Ms::Vector{MPS};  weight::Float64, recalcEnv::Bool = true)
updateH!(sysenv::StateEnvs{ProjMPOSum_MPS}, H::Vector{MPO}, Ms::Vector{MPS}; weight::Float64, recalcEnv::Bool = true)
updateH!(sysenv::StateEnvs{ProjCouplingModel_MPS}, H::CouplingModel, Ms::Vector{MPS}; weight::Float64, recalcEnv::Bool = true)
nsite(sysenv::StateEnvs)
set_nsite!(sysenv::StateEnvs, nsite::Int)
position!(sysenv::StateEnvs, pos::Int)
TenNetLib.product(sysenv::StateEnvs, v::ITensor)
Base.copy(sysenv::StateEnvs)
Base.length(sysenv::StateEnvs)
```