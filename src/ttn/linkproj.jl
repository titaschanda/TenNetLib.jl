
#################################################################################
"""
    mutable struct LinkProjTTN{T}
        tensors::Dict{LinkTypeTTN, ITensor}
        M::TTN{T}
    end

Holds projectors at different links of the TTN. Required to calculate excited state or
overlap with another state.

 - `tensors::Dict{LinkTypeTTN, ITensor}`: Tensors are each link.
 - `M::TTN`: The TTN to project.
"""
mutable struct LinkProjTTN{T}
    tensors::Dict{LinkTypeTTN, ITensor}    
    M::TTN{T}
end
    
#################################################################################

"""
    Base.getindex(env::LinkProjTTN, node1::Int2, node2::Int2)
    Base.getindex(env::LinkProjTTN, link::LinkTypeTTN)

Returns the projector tensor at a link.
"""
Base.getindex(env::LinkProjTTN, node1::Int2, node2::Int2) =
    Base.getindex(env.tensors, LinkTypeTTN(node1, node2))

Base.getindex(env::LinkProjTTN, link::LinkTypeTTN) =
    Base.getindex(env.tensors, link)

#################################################################################

"""
    Base.setindex!(env::LinkProjTTN, tensors::ITensor, node1::Int2, node2::Int2)
    Base.setindex!(env::LinkProjTTN, tensors::ITensor, link::LinkTypeTTN)

Sets the projector tensor at a link.
"""
Base.setindex!(env::LinkProjTTN, tensors::ITensor, node1::Int2, node2::Int2) =
    Base.setindex!(env.tensors, tensors, LinkTypeTTN(node1, node2))

Base.setindex!(env::LinkProjTTN, tensors::ITensor, link::LinkTypeTTN) =
    Base.setindex!(env.tensors, tensors, link)

#################################################################################

"""
    function move_linkproj_to_next!(env::LinkProjTTN, psi::TTN,
                                    node::Int2, nextnode::Int2)::Nothing

Moves the projectors to a neighboring node.
"""
function move_linkproj_to_next!(env::LinkProjTTN, psi::TTN,
                                node::Int2, nextnode::Int2)::Nothing
    
    if !isvalidnode(psi, node) || !isvalidnode(psi, nextnode)
        error("`move_linkproj_to_next!()`: Input position is ill-defined !!")
    end
    if !isneighbor(psi, node, nextnode)
        error("`move_linkproj_to_next!()`: Input nodes are not neighbors !!")
    end

    next_link, prev_links  = _get_links(psi, node, nextnode)    
    tens = ITensor[env[x] for x in prev_links if haskey(env.tensors, x)]
    
    phiA = psi[node]    
    env[next_link] = dag(prime(env.M[node]; tags = "Link"))

    if length(tens) > 0
        env[next_link] = ITensors.contract(env[next_link],
                                           tens...;
                                           sequence = ITensors.default_sequence())
    end
    env[next_link] *= phiA

    return nothing
end

#################################################################################

"""
    function move_linkproj!(env::LinkProjTTN, psi::TTN,
                            source_node::Int2,
                            destination_node::Int2;
                            kwargs...)::Nothing

Moves the projectors from `source_node` to `destination_node`.
"""
function move_linkproj!(env::LinkProjTTN, psi::TTN,
                        source_node::Int2,
                        destination_node::Int2;
                        kwargs...)::Nothing

    if !isvalidnode(psi, source_node) || !isvalidnode(psi, destination_node)
        error("`move_linkproj!()`: Input position is ill-defined !!")
    end

    if source_node == destination_node
        return nothing
    end

    to_omit::Union{Int2, Nothing} = get(kwargs, :node_to_skip, nothing)    
    path = shortest_path(psi.graph, source_node, destination_node)    
    for dummy = 1 : length(path) - 1
        node1 = path[dummy]
        node2 = path[dummy+1]        
        if node1 == to_omit
            continue
        end
        move_linkproj_to_next!(env, psi, node1, node2)
    end
    return nothing
end

#################################################################################

"""
    function Base.eltype(env::LinkProjTTN, psi::TTN)

Element type (i.e., `Float64` or `ComplexF64`).
"""
function Base.eltype(env::LinkProjTTN, psi::TTN)
    oc = orthocenter(psi)
    links = _get_links(psi, oc)
    elt = Bool
    for ll in links
        if haskey(env.tensors, ll)
            elt = promote_type(elt, eltype(env[ll]))
        end
    end
    return elt
end
    
#################################################################################

"""
    function contract(env::LinkProjTTN, psi::TTN, v::ITensor)::ITensor

Returns the contraction of the projector with input `v`. The contraction is performed at
the orthogonality center of `psi`.
"""
function contract(env::LinkProjTTN, psi::TTN, v::ITensor)::ITensor
    
    oc = orthocenter(psi)
    links = _get_links(psi, oc)
    tens = ITensor[env[x] for x in links if haskey(env.tensors, x)]

    Hv = ITensors.contract(dag(prime(env.M[oc]; tags = "Link")),
                           tens...;
                           sequence = ITensors.default_sequence())
    Hv *= v
    return Hv
end

#################################################################################

"""
    function product(env::LinkProjTTN, psi::TTN, v::ITensor)::ITensor

Returns the product of the projector with input `v`. The contraction is performed at
the orthogonality center of `psi`.
"""
function product(env::LinkProjTTN, psi::TTN, v::ITensor)::ITensor
    Pv = contract(env, psi, v) * dag(contract(env, psi, ITensor(true)))
    if order(Pv) != order(v)
        error(
            string(
                "The order of the LinkProjTTN-ITensor product P*v ", 
                "is not equal to the order of the ITensor v, ",
                "this is probably due to an index mismatch."
            ),
        )
    end
    
    return noprime(Pv)
end

#################################################################################

"""
    function LinkProjTTN(psi::TTN{T}, M::TTN{T})::LinkProjTTN where T

Constructor for the `LinkProjTTN`.
"""
function LinkProjTTN(psi::TTN{T}, M::TTN{T})::LinkProjTTN where T

    @assert siteinds(psi) == siteinds(M)
    @assert psi.graph == M.graph
        
    if !isvalidnode(psi, orthocenter(psi))
        error("`LinkProjTTN()`: TTN does not have a proper orthogonality center !!")
    end

    if hasqns(psi)
        qnnode_psi = find_qnnode(psi)
        qnind_psi = only(filter(hastags("QN"), inds(psi[qnnode_psi])))

        qnnode_M = find_qnnode(M)
        qnind_M = only(filter(hastags("QN"), inds(M[qnnode_M])))

        if (space(qnind_psi) != space(qnind_M))
            error("`LinkProjTTN()`: TTNs do have same global QN !!")
        end
        
        replaceind!(psi[qnnode_M], qnind_psi, qnind_M)
    end
        
    env = LinkProjTTN(Dict{LinkTypeTTN, ITensor}(), M)
    
    oc = orthocenter(psi)
    path = nodes_from_bfs(psi.graph, oc; reverse=true)
    for dummy = 1 : length(path) - 1
        node1 = path[dummy]
        node2 = nextnode_in_path(psi.graph, node1, oc)
        move_linkproj_to_next!(env, psi, node1, node2)
    end

    return env
end

#################################################################################
