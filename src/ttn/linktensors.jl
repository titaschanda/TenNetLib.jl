
#################################################################################

"""
    const LinkTensorsTTN = Dict{LinkTypeTTN, IDTensors}

Holds the environment tensors at different links of the TTN.
"""
const LinkTensorsTTN = Dict{LinkTypeTTN, IDTensors}

#################################################################################

"""
    Base.getindex(env::LinkTensorsTTN, node1::Int2, node2::Int2)

Returns the environment tensors at a link.
"""
Base.getindex(env::LinkTensorsTTN, node1::Int2, node2::Int2) =
    Base.getindex(env, LinkTypeTTN(node1, node2))

#################################################################################

"""
    Base.setindex!(env::LinkTensorsTTN, tensors::IDTensors, node1::Int2, node2::Int2)

Sets the environment tensors at a link.
"""
Base.setindex!(env::LinkTensorsTTN, tensors::IDTensors, node1::Int2, node2::Int2) =
    Base.setindex!(env, tensors, LinkTypeTTN(node1, node2))

#################################################################################

function _collect_link_tensors(env::LinkTensorsTTN,
                               links::Vector{LinkTypeTTN},
                               hasqn::Bool
                               )::Tuple{Dict{IDType, Vector{ITensor}}, Bool}    
    removeqn = false
    idtens = Dict{IDType, Vector{ITensor}}()
    for link in links
        if haskey(env, link)
            for term in env[link]
                !hasqns(term.second) && hasqn && (removeqn = true)
                if !haskey(idtens, term.first)
                    idtens[term.first] = ITensor[term.second]
                else
                    push!(idtens[term.first], term.second)
                end
            end
        end
    end

    return idtens, removeqn
end

#################################################################################

"""
    function move_linktensors_to_next!(env::LinkTensorsTTN, psi::TTN,
                                       node::Int2, nextnode::Int2)::Nothing

Moves the environment tensors to a neighboring node.
"""
function move_linktensors_to_next!(env::LinkTensorsTTN, psi::TTN,
                                   node::Int2, nextnode::Int2)::Nothing
    
    if !isvalidnode(psi, node) || !isvalidnode(psi, nextnode)
        error("`move_linktensors_to_next!()`: Input position is ill-defined !!")
    end
    if !isneighbor(psi, node, nextnode)
        error("`move_linktensors_to_next!()`: Input nodes are not neighbors !!")
    end

    next_link, prev_links  = _get_links(psi, node, nextnode)
    idtens, removeqn = _collect_link_tensors(env, prev_links, hasqns(psi))
        
    length(idtens) == 0 && (return nothing)
    conditional_removeqns(A) = removeqn && hasqns(A) ? removeqns(A) : A

    env[next_link] = IDTensors()    
    nextind = conditional_removeqns(commonind(psi[node], psi[nextnode]))
    local_tensor = ITensor()
    idkeys = collect(keys(idtens))

    using_threaded_loop() && (mutex = Threads.SpinLock())
    @threaded_loop for id in idkeys
        phi = conditional_removeqns(psi[node])
        inds_to_prime = typeof(nextind)[commonind(phi, x) for x in idtens[id]]
        push!(inds_to_prime, nextind)
        phidag = conditional_removeqns(dag(prime(phi, inds_to_prime)))

        phidag = ITensors.contract(phidag,
                                   conditional_removeqns.(idtens[id])...;
                                   sequence = ITensors.default_sequence())        
        phidag *= phi

        if using_threaded_loop()
            lock(mutex) do
                if order(phidag) > 2
                    env[next_link][id] = phidag                
                else
                    local_tensor += phidag
                end  
            end
        else
            if order(phidag) > 2
                env[next_link][id] = phidag                
            else
                local_tensor += phidag
            end 
        end
    end

    if order(local_tensor) != 0
        newid = gen_rand_id()
        env[next_link][newid] = local_tensor
    end
    return nothing
end

#################################################################################

"""
    function move_linktensors!(env::LinkTensorsTTN, psi::TTN,
                               source_node::Int2,
                               destination_node::Int2;
                               kwargs...)::Nothing

Moves the environment tensors from `source_node` to `destination_node`.
"""
function move_linktensors!(env::LinkTensorsTTN, psi::TTN,
                           source_node::Int2,
                           destination_node::Int2;
                           kwargs...)::Nothing

    if !isvalidnode(psi, source_node) || !isvalidnode(psi, destination_node)
        error("`move_linktensors!()`: Input position is ill-defined !!")
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
        move_linktensors_to_next!(env, psi, node1, node2)
    end
    return nothing
end

#################################################################################

"""
    function Base.eltype(env::LinkTensorsTTN, psi::TTN)

Element type (i.e., `Float64` or `ComplexF64`).
"""
function Base.eltype(env::LinkTensorsTTN, psi::TTN)
    oc = orthocenter(psi)
    links = _get_links(psi, oc)
    elt = Bool
    for ll in links
        if haskey(env, ll)
            elt = promote_type(elt, eltype(env[ll]))
        end
    end
    return elt
end
    
#################################################################################

"""
    function product(env::LinkTensorsTTN, psi::TTN, v::ITensor)::ITensor

Returns the product of the environment tensors with input `v`. The contraction is performed at
the orthogonality center of `psi`.
"""
function product(env::LinkTensorsTTN, psi::TTN, v::ITensor)::ITensor   
    
    oc = orthocenter(psi)
    links = _get_links(psi, oc)

    idtens, removeqn = _collect_link_tensors(env, links, hasqns(v))    
    conditional_removeqns(A) = removeqn && hasqns(A) ? removeqns(A) : A
    
    sum_tensor = ITensor()    
    idkeys = collect(keys(idtens))

    using_threaded_loop() && (mutex = Threads.SpinLock())   
    @threaded_loop for key in idkeys
        
        Hv = ITensors.contract(conditional_removeqns(v),
                               conditional_removeqns.(idtens[key])...;
                               sequence = ITensors.default_sequence())

        if using_threaded_loop()
            lock(mutex) do
                sum_tensor += noprime(Hv)
            end
        else
            sum_tensor += noprime(Hv)
        end
    end
    
    if order(sum_tensor) != order(v)
        error(
            string(
                "The order of the LinkTensorsTTN-ITensor product P*v ", 
                "is not equal to the order of the ITensor v, ",
                "this is probably due to an index mismatch."
            ),
        )
    end
            
    return sum_tensor
end

#################################################################################

"""
    function LinkTensorsTTN(psi::TTN, M::CouplingModel)::LinkTensorsTTN
    function LinkTensorsTTN(psi::TTN, optens::Vector{ITensor})::LinkTensorsTTN
    function LinkTensorsTTN(psi::TTN, oppairs::Vector{Pair{String, Int}};
                            isfermions::Bool = true)::LinkTensorsTTN

Constructors for the `LinkTensorsTTN`.
"""
function LinkTensorsTTN(psi::TTN, M::CouplingModel)::LinkTensorsTTN

    @assert siteinds(psi) == siteinds(M)
    
    if !isvalidnode(psi, orthocenter(psi))
        error("`LinkTensorsTTN()`: TTN does not have a proper orthogonality center !!")
    end

    env = LinkTensorsTTN()
    
    nodelist = Set{Int2}()
    for n in 1 : numsites(psi)
        node = find_sitenode(psi, n)

        if length(M[n]) > 0
            env[(0, n), node] = M[n]
            push!(nodelist, node)
        end
    end

    oc = orthocenter(psi)
    path = nodes_from_bfs(psi.graph, oc, nodelist; reverse=true)
    for dummy = 1 : length(path) - 1
        node1 = path[dummy]
        node2 = nextnode_in_path(psi.graph, node1, oc)
        move_linktensors_to_next!(env, psi, node1, node2)
    end

    return env
end

#################################################################################

function LinkTensorsTTN(psi::TTN, optens::Vector{ITensor})::LinkTensorsTTN

    if !isvalidnode(psi, orthocenter(psi))
        error("`LinkTensorsTTN()`: TTN does not have a proper orthogonality center !!")
    end

    #************************************************************
    
    oppairs = Vector{Pair{ITensor, Int}}()
    removeqn = false
    for o in optens
        pos = findsite(psi, o)
        pass = hasind(o, siteind(psi, pos)) && hasind(o, siteind(psi, pos)') &&
            ndims(o) == 2
        if !pass
            error("`LinkTensorsTTN()`: Error in operator tensors !!")
        end
        !hasqns(o) && hasqns(psi) && (removeqn = true)
        push!(oppairs, o => pos)
    end
    sort!(oppairs, by = x -> x.second)
    conditional_removeqns(A) = removeqn && hasqns(A) ? removeqns(A) : A
    
    #************************************************************

    opdict = Dict{Int, ITensor}()
    for o in oppairs
        if !haskey(opdict, o.second)
            opdict[o.second] = conditional_removeqns(o.first)
        else
            opdict[o.second] = prime(opdict[o.second])
            opdict[o.second] *= conditional_removeqns(o.first)
            opdict[o.second] = mapprime(opdict[o.second], 2, 1)
        end
    end

    positions = sort(collect(keys(opdict)))
    for ii = 2 : length(positions)
        dummyindex = !removeqn ? 
            Index(QN() => 1; tags = "OpLink") : 
            Index(1; tags = "OpLink")
        opdict[positions[ii-1]] *= onehot(dummyindex => 1)
        opdict[positions[ii]] *= onehot(dag(dummyindex) => 1)
    end
    
    #************************************************************
    
    env = LinkTensorsTTN()
    id  = gen_rand_id()
    
    for n in keys(opdict)
        node = find_sitenode(psi, n)        
        env[(0, n), node] = IDTensors(id => opdict[n])
    end

    #************************************************************

    oc = orthocenter(psi)
    path = nodes_from_bfs(psi.graph, oc,
                          find_sitenodes(psi,
                                         collect(keys(opdict)));
                          reverse=true)
    for dummy = 1 : length(path) - 1
        node1 = path[dummy]
        node2 = nextnode_in_path(psi.graph, node1, oc)
        move_linktensors_to_next!(env, psi, node1, node2)
    end

    return env
end

#################################################################################

function LinkTensorsTTN(psi::TTN, oppairs::Vector{Pair{String, Int}};
                        isfermions::Bool = true)::LinkTensorsTTN

    if isfermions
        perm, oppairs = bosonize(oppairs, siteinds(psi))
    end
    optens = ITensor[op(x.first, siteind(psi, x.second)) for x in oppairs]
    isfermions && (optens[begin] *= perm)
    return LinkTensorsTTN(psi, optens)
end

#################################################################################

