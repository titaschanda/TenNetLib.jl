# Perform local updates

At the lowest-level of abstraction, TenNetLib.jl allows for updating the `StateEnvsTTN` for each sites manually.

Skip this part if you want to avoid lower-level abstraction.

```@docs
update_position!(sysenv::StateEnvsTTN, solver, node::Int2; time_step::Union{Float64, ComplexF64, Nothing}, normalize::Bool, maxdim::Int, mindim::Int, cutoff::Float64, svd_alg::String, kwargs...)
```

TenNetLib.jl implements the subspace expansion method described in
[SciPost Phys. Lect. Notes 8 (2019)] (https://scipost.org/10.21468/SciPostPhysLectNotes.8) to increase the bond dimension between two neighboring nodes.

```@docs
subspace_expand!(psi::TTN, node::Int2, nextnode::Int2, max_expand_dim::Int, noise::Float64)
```