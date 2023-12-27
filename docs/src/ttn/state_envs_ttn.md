# `StateEnvsTTN`: A container to store the TTN and its environments

Similar to the setup for MPS, tt the lowest-level of abstraction, TenNetLib.jl defines
`StateEnvs` to hold a TTN and its environments to be modified in place.

Skip this part if you want to avoid lower-level abstraction.

```@docs
StateEnvsTTN
StateEnvsTTN(psi::TTN, M::CouplingModel)
StateEnvsTTN(psi::TTN,M::CouplingModel,psis::Vector{TTN};  weight::Float64)
getpsi(sysenv::StateEnvsTTN)
getenv(sysenv::StateEnvsTTN)
position!(sysenv::StateEnvsTTN, node::Int2; maxdim::Int = typemax(Int), mindim::Int = 1, cutoff::Float64 = Float64_threshold(), svd_alg::String = "divide_and_conquer", normalize::Bool = false, node_to_skip::Union{Int2, Nothing} = nothing)
TenNetLib.product(sysenv::StateEnvsTTN, v::ITensor)
Base.copy(sysenv::StateEnvsTTN)
```