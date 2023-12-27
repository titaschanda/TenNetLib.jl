# `CouplingModel`

TenNetLib.jl degines as struct, called the `CouplingModel`, to store the Hamiltonian terms. 
In case of MPS based algorithms, `CouplingModel` can replace `MPO` without modifying
rest of the code. For Tree Tensor Network (TTN) codes, only `CouplingModel` can be used.
Different elements of `CouplingModel` are contracted in parallel.

```@docs
CouplingModel
CouplingModel(os::OpStrings{T1}, sites::Vector{Index{T2}}; merge::Bool = true,maxdim::Int = typemax(Int), mindim::Int = 1, cutoff::Float64 = Float64_threashold(),svd_alg::String = "divide_and_conquer", chunksize::Int = 12) where {T1 <: Number, T2}
CouplingModel(os::OpStrings{T1}, mpo::MPO; merge::Bool = true,maxdim::Int = typemax(Int), mindim::Int = 1, cutoff::Float64 = Float64_threashold(),svd_alg::String = "divide_and_conquer", chunksize::Int = 12) where {T1 <: Number}
CouplingModel(mpos::MPO...)
Base.length(model::CouplingModel)
Base.copy(model::CouplingModel)
Base.getindex(model::CouplingModel, n)
ITensors.siteinds(model::CouplingModel)
```