# The `TTN` object

```@docs
TTN
ITensorMPS.siteinds(ttn::TTN)
ITensorMPS.siteind(ttn::TTN, n::Int)
numsites(ttn::TTN)
getgraph(ttn::TTN)
Base.getindex(ttn::TTN, node::Int2)
Base.setindex!(ttn::TTN, tensor::ITensor, node::Int2)
orthocenter(ttn::TTN)
Base.copy(ttn::TTN)
isvalidnode(ttn::TTN, node::Int2)
isneighbor(ttn::TTN, node1::Int2, node2::Int2)
ITensorMPS.findsite(ttn::TTN, is)
ITensorMPS.findsites(ttn::TTN, is)
findnode(ttn::TTN, is)
findnodes(ttn::TTN, is)
find_sitenode(ttn::TTN, n::Int)
find_sitenodes(ttn::TTN, ns)
ITensorMPS.maxlinkdim(ttn::TTN)
LinearAlgebra.normalize!(ttn::TTN)
```

!!! info
    To conserve global symmetry in TTN, TenNetLib.jl uses a dummy index of dimension one whose flux
    fixes the superselection sector. This dummy index is attached to one of the tensors of TTN,
    usually the eccentric central node. Unlike MPS, the fluxes of all tensors in TTN are zero.

```@docs
ITensors.hasqns(ttn::TTN)
find_qnnode(ttn::TTN{ITensors.QNBlocks})
moveisometry_to_next!(ttn::TTN, node1::Int2, node2::Int2; kwargs...)
isometrize_full!(ttn::TTN, node::Int2; kwargs...)
isometrize!(ttn::TTN, node::Int2; kwargs...)
```