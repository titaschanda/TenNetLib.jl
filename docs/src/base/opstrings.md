# `OpStrings`

TeNLib provides an alternative, called `OpStrings`, to ITensors' `OpSum` to automatically
construct Hamiltonains / operators. TeNLib's own `CouplingModel` is built from `OpStrings` and
is not compatible with `OpSum`.

## `OpString`

`OpStrings` is basically a vector of `OpString` objects (notice the difference in 's' at the end)

```@docs
OpString
TeNLib.coefficient(opstr::OpString)
operators(opstr::OpString)
minsite(opstr::OpString)
maxsite(opstr::OpString)
removeIds(opstr::OpString{T}) where {T <: Number}
bosonize(opstr::OpString{T1}, sites::Vector{Index{T2}}) where {T1 <: Number, T2}
```

## `OpStrings`

```@docs
OpStrings{T}
removeIdsZeros(os::OpStrings{T}) where {T <: Number}
bosonize(os::OpStrings{T1}, sites::Vector{Index{T2}}) where {T1 <: Number, T2}
mergeterms(os::OpStrings{T}) where T <: Number
```

## `MPO` from `OpStrings`
```@docs
ITensors.MPO(os::OpStrings{T1}, sites::Vector{Index{T2}}; maxdim::Int = typemax(Int), mindim::Int = 1, cutoff::Float64 = Float64_threashold(), svd_alg::String = "divide_and_conquer", chunksize::Int = 12) where {T1 <: Number, T2}
```