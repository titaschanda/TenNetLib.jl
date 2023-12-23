
#################################################################################

"""
    mutable struct TTN{T}
        sites::Vector{Index{T}}
        graph::Graph{Int2}
        tensors::Dict{Int2, ITensor}
        orthocenter::Int2
    end

Tree Tensor Network (TTN) object.

Each tensor is indexed by the tuple of integers (`::Int2`), `(ll, nn)`,
where (usually) `ll` denotes the layer index and `nn` denotes the tensor index at each layer.

 - `sites::Vector{Index}`: Physical site `Index`s.
 - `graph::Graph{Int2}`: Underlying structure of the network defined by `Graph{Int2}`.
 - `tensors::Dict{Int2, ITensor}`: Tensors in the TTN.
 - `orthocenter::Int2`: Orthogonality center of the TTN.
"""
mutable struct TTN{T}
    sites::Vector{Index{T}}
    graph::Graph{Int2}
    tensors::Dict{Int2, ITensor}
    orthocenter::Int2
end

#################################################################################

"""
    function ITensors.siteinds(ttn::TTN)

Returns the site `Index`s of the `TTN`.
"""
ITensors.siteinds(ttn::TTN) = ttn.sites

#################################################################################

"""
    function ITensors.siteind(ttn::TTN, n::Int)

Returns the `n`th site `Index` of the `TTN`.
"""
ITensors.siteind(ttn::TTN, n::Int) = ttn.sites[n]

#################################################################################

"""
    function numsites(ttn::TTN)

Returns the number of physical sites in the `TTN`.
"""
numsites(ttn::TTN) = length(ttn.sites)

#################################################################################

"""
    function getgraph(ttn::TTN)

Returns (shallow copy of) the underlying graph of the `TTN`.
"""
getgraph(ttn::TTN) = copy(ttn.graph)

#################################################################################

"""
    function Base.getindex(ttn::TTN, node::Int2)
    function Base.getindex(ttn::TTN, ll::Int, nn::Int)

Returns the tensor at the `node::Int2` = `(ll, nn)`.
"""
Base.getindex(ttn::TTN, node::Int2) = Base.getindex(ttn.tensors, node)
Base.getindex(ttn::TTN, ll::Int, nn::Int) = Base.getindex(ttn, (ll, nn))


#################################################################################

"""
    function Base.setindex!(ttn::TTN, tensor::ITensor, node::Int2)
    function Base.setindex!(ttn::TTN, tensor::ITensor, ll::Int, nn::Int)

Sets the tensor at the `node::Int2` = `(ll, nn)`. If the `node` is not the orthogonality center,
orthogonality center is set to `(typemin(Int), typemin(Int))`.
"""
function Base.setindex!(ttn::TTN, tensor::ITensor, node::Int2)
    Base.setindex!(ttn.tensors, tensor, node)
    if orthocenter(ttn) != node
        ttn.orthocenter = Int2()
    end
end

Base.setindex!(ttn::TTN, tensor::ITensor, ll::Int, nn::Int) =
    Base.setindex!(ttn, tensor, (ll, nn))

#################################################################################

"""
    function orthocenter(ttn::TTN)

Returns the orthogonality center of the `TTN`.
"""
orthocenter(ttn::TTN) = ttn.orthocenter

#################################################################################

"""
    Base.copy(ttn::TTN)

Shallow copy of `TTN`.
"""
Base.copy(ttn::TTN) = TTN(Base.copy(ttn.sites),
                          Base.copy(ttn.graph),
                          Base.copy(ttn.tensors),
                          ttn.orthocenter)

#################################################################################

"""
    function isvalidnode(ttn::TTN, node::Int2)

Checks whether the input `node::Int2` is a valid node in the `TTN` graph.
"""
isvalidnode(ttn::TTN, node::Int2) = node in ttn.graph.nodes

#################################################################################

"""
    function isneighbor(ttn::TTN, node1::Int2, node2::Int2)

Checks whether `node1` and `node2` are neighboring nodes in the TTN.
"""
isneighbor(ttn::TTN, node1::Int2, node2::Int2)::Bool =
    isneighbor(ttn.graph, node1, node2)

#################################################################################


"""
    function ITensors.findsite(ttn::TTN, is)

Returns the first site of the TTN that has at least one `Index` in
common with the `Index` or `ITensor` or their collection `is`.
"""
ITensors.findsite(ttn::TTN, is) = findfirst(hascommoninds(is), ttn.sites)


"""
    function ITensors.findsites(ttn::TTN, is)

Returns the sites of the TTN that have `Index`s in
common with the `Index` or `ITensor` or their collection `is`.
"""
ITensors.findsites(ttn::TTN, is) = findall(hascommoninds(is), ttn.sites)

#################################################################################

"""
    function findnode(ttn::TTN, is)

Returns the first node of the TTN that has at least one `Index` in
common with the `Index` or `ITensor` or their collection `is`.
"""
findnode(ttn::TTN, is) = findfirst(hascommoninds(is), ttn.tensors)


"""
    function findnodes(ttn::TTN, is)

Returns the nodes of the TTN that have `Index`s in
common with the `Index` or `ITensor` or their collection `is`.
"""
findnodes(ttn::TTN, is) = findall(hascommoninds(is), ttn.tensors)

#################################################################################

"""
    function find_sitenode(ttn::TTN, n::Int)

Return node of the TTN that is associated with the site `n`. 
"""
find_sitenode(ttn::TTN, n::Int) = findfirst(hascommoninds(ttn.sites[n]), ttn.tensors)


"""
    function find_sitenodes(ttn::TTN, ns)

Returns nodes of the TTN that are associated with the sites `ns` (a collection of `Int`s). 
"""
find_sitenodes(ttn::TTN, ns) = map(x -> find_sitenode(ttn, x), ns)

#################################################################################

"""
    function ITensors.maxlinkdim(ttn::TTN)

Returns the maximum bond/link dimension of the TTN.
"""
function ITensors.maxlinkdim(ttn::TTN)
    md = 1
    for node in ttn.graph.nodes
        l = filter(hastags("Link"), inds(ttn[node]))
        linkdim = length(l) == 0 ? 1 : maximum(Int[dim(x) for x in l])
        md = max(md, linkdim)
    end
    return md
end

#################################################################################

"""
    function LinearAlgebra.normalize!(ttn::TTN)

Normalizes the TTN. The TTN must have well-defined orthogonality center.
"""
function LinearAlgebra.normalize!(ttn::TTN)::TTN
    if !isvalidnode(ttn, orthocenter(ttn))
        error("`normalize!()`: TTN does not have a proper orthogonality center !!")
    end
    normalize!(ttn[orthocenter(ttn)])
    return ttn
end

#################################################################################

"""
    function ITensors.hasqns(ttn::TTN)

Checks whether the TTN has QNs.
"""
ITensors.hasqns(ttn::TTN) = hasqns(ttn.sites)

#################################################################################

"""
    function find_qnnode(ttn::TTN)

Returns the node where the a dummy QN Index is attached to fix the QN sector.
Returns `nothing` if the TTN does not have QN.
"""
find_qnnode(ttn::TTN{QNBlocks}) =
    findfirst(hastags("QN"), ttn.tensors)

find_qnnode(ttn::TTN{Int}) = nothing


#################################################################################

"""
    function moveisometry_to_next!(ttn::TTN, node1::Int2, node2::Int2; kwargs...)

Moves the isometry / orthogonality center of the TTN from `node1::Int2` to the
neighboring `node2::Int`. `node1` and `node2` must be neighboring nodes.
`node1` must be the orthogonality center unless `ignore_orthocenter = true`.

**Note**: Does not normalize the TTN afterwards.

#### Named arguments and their default values.
 - `ignore_orthocenter::Bool` = false: If set to `true`, `node1` is not required to
   be the orthogonality center.
 - `maxdim::Int = typemax(Int)`: Maximum bond dimension after SVD truncation.
 - `mindim::Int = 1`: Minimum bond dimension after SVD truncation.
 - `cutoff::Float64 = Float64_threshold()`: Cutoff for SVD truncation.
 - `svd_alg::String = "divide_and_conquer"`.
"""
function moveisometry_to_next!(ttn::TTN, node1::Int2, node2::Int2; kwargs...)
    
    ignoreOC::Bool = get(kwargs, :ignore_orthocenter, false)
    
    if !isvalidnode(ttn, node1) || !isvalidnode(ttn, node2)
        error("`moveisometry_to_next!()`: Input position is ill-defined !!")
    end
    if orthocenter(ttn) != node1 && !ignoreOC
        error(string("`moveisometry_to_next!()`: `orthocenter` does not match with the ",
                     "input `node1 = $node1 !!"))
    end
    if !isneighbor(ttn, node1, node2)
        error("`moveisometry_to_next!()`: Input nodes are not neighbors !!")
    end
    
    maxdim::Int = get(kwargs, :maxdim, typemax(Int))
    mindim = get(kwargs, :mindim, 1)
    cutoff::Float64 = get(kwargs, :cutoff, Float64_threshold())
    svd_alg::String = get(kwargs, :svd_alg, "divide_and_conquer")
    
    A = ttn[node1]
    B = ttn[node2]
    comind = commonind(A, B)
    uinds = uniqueinds(A, B)
    tag = tags(comind)

    if cutoff > 0.0 && dim(comind) > maxdim
        U, S, V = svd(A, uinds;
                      maxdim = maxdim,
                      mindim = mindim,
                      cutoff = cutoff,
                      alg = svd_alg,
                      lefttags = tag)
        ttn[node1] = U
        ttn[node2] *= S*V
    else
        U, V = qr(A, uinds,
                  tags = tag)
        ttn[node1] = U
        ttn[node2] *= V
    end
    
    ttn.orthocenter = node2
    return nothing
end

#################################################################################

"""
    function isometrize_full!(ttn::TTN, node::Int2; kwargs...)

Isometrizes the TTN from scratch with respect to the orthogonality center `node::Int2`.

#### Named arguments and their default values.
 - `normalize::Bool = true`: Whether to normalize the TTN afterwards.
 - `maxdim::Int = typemax(Int)`: Maximum bond dimension after SVD truncation.
 - `mindim::Int = 1`: Minimum bond dimension after SVD truncation.
 - `cutoff::Float64 = Float64_threshold()`: Cutoff for SVD truncation.
 - `svd_alg::String = "divide_and_conquer"`. 
"""
function isometrize_full!(ttn::TTN, node::Int2; kwargs...)
    if !isvalidnode(ttn, node)
        error("`isometrize_full!()`: Input position is ill-defined !!")
    end

    bfs_path = nodes_from_bfs(ttn.graph, node; reverse=true)
      
    for ii = 1 : length(bfs_path)-1
        nextnode = nextnode_in_path(ttn.graph, bfs_path[ii], node)
        moveisometry_to_next!(ttn, bfs_path[ii], nextnode; ignore_orthocenter = true, kwargs...)
    end
    
    normalize::Bool = get(kwargs, :normalize, true)
    normalize && normalize!(ttn)

    return nothing
end

#################################################################################

"""
    function isometrize!(ttn::TTN, node::Int2; kwargs...)

Moves the isometry / orthogonality center of the TTN to a new center `node`.
If the TTN does not have proper orthoginality center, Isometrizes the TTN from scratch.

#### Named arguments and their default values.
 - `normalize::Bool = true`: Whether to normalize the TTN afterwards.
 - `maxdim::Int = typemax(Int)`: Maximum bond dimension after SVD truncation.
 - `mindim::Int = 1`: Minimum bond dimension after SVD truncation.
 - `cutoff::Float64 = Float64_threshold()`: Cutoff for SVD truncation.
 - `svd_alg::String = "divide_and_conquer"`. 
"""
function isometrize!(ttn::TTN, node::Int2; kwargs...)
    if !isvalidnode(ttn, node)
        error("`isometrize!()`: Input position is ill-defined !!")
    end
    
    if !isvalidnode(ttn, orthocenter(ttn))
        isometrize_full(ttn, node; kwargs...)
    end
  
    normalize::Bool = get(kwargs, :normalize, true)
    
    if node == orthocenter(ttn)
        normalize && normalize!(ttn)
        return nothing
    
    else
        direc = shortest_path(ttn.graph, orthocenter(ttn), node)
        for ii = 1 : length(direc) - 1
            moveisometry_to_next!(ttn, direc[ii], direc[ii+1]; kwargs...)
        end
        normalize && normalize!(ttn)
        return nothing
    end
end

#################################################################################
